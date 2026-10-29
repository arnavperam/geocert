"""Run one frozen scientific stage in a fresh process and record its resources."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import resource
import runpy
import sys
import time
import traceback


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--script', type=Path, required=True)
    parser.add_argument('--record', type=Path, required=True)
    args = parser.parse_args()
    started = time.time()
    timer = time.perf_counter()
    status, error, code = 'passed', None, 0
    sys.argv = [str(args.script)]
    try:
        runpy.run_path(str(args.script), run_name='__main__')
    except SystemExit as exc:
        code = int(exc.code or 0)
        if code:
            status, error = 'failed', f'SystemExit({code})'
    except BaseException as exc:
        status, error, code = 'failed', f'{type(exc).__name__}: {exc}', 1
        traceback.print_exc()
    usage = resource.getrusage(resource.RUSAGE_SELF)
    record = {
        'script': args.script.name,
        'status': status,
        'error': error,
        'exit_code': code,
        'started_unix': started,
        'ended_unix': time.time(),
        'wall_seconds': time.perf_counter() - timer,
        'user_cpu_seconds': usage.ru_utime,
        'system_cpu_seconds': usage.ru_stime,
        'peak_process_rss_kib': usage.ru_maxrss,
        'rss_semantics': 'Linux ru_maxrss, this fresh process including imports and output serialization',
        'python': sys.version,
        'platform': platform.platform(),
        'cpu_affinity': sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else None,
        'thread_environment': {k: os.environ.get(k) for k in
                               ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']},
    }
    args.record.write_text(json.dumps(record, indent=2) + '\n')
    return code


if __name__ == '__main__':
    raise SystemExit(main())
