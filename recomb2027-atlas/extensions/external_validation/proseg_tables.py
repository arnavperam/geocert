"""Build the complete setting inventory and supplementary tables after all runs."""
from pathlib import Path
import argparse
import csv
import json


def build(analysis, output, tex):
    rows = []
    background = []
    for donor, prefix, number in [('LUNG', 'roi', 2), ('IDC', 'idc', 3), ('ILC', 'ilc', 3)]:
        for index in range(1, number + 1):
            region = analysis / donor / f'{prefix}_{index:02d}'
            comparison = json.loads((region / 'results/proseg_comparison.json').read_text())
            records = {r['run']: r for r in comparison['runs']}
            bg = {'core': f'{donor} {index}'}
            for side in [200, 800]:
                for variant in ['default', 'no_factorization']:
                    name = f'{side}_{variant}'
                    directory = region / 'proseg' / name
                    run = json.loads((directory / 'run.json').read_text())
                    row = {'donor': donor, 'region': region.name, 'input_um': side,
                           'factorization': variant, 'returncode': run['returncode'], 'seconds': run['seconds']}
                    if run['returncode'] == 0:
                        audit = json.loads((directory / 'audit.json').read_text())
                        assert audit['failures'] == 0
                        assert (directory / 'core_assignments.json.gz').exists()
                        record = records[name]
                        row.update(record)
                        row['background_percent'] = 100 * record['background_core_RNA'] / record['common_retained_core_RNA']
                        row['reach_percent'] = 100 * record['changed_counts_reachable_at4'] / record['changed_counts']
                        row['background_related_percent'] = 100 * record['changed_counts_involving_background'] / record['changed_counts']
                        row['native_count_matrix_audit_failures'] = audit['failures']
                        row['optional_csv_export_complete'] = audit['optional_csv_export_complete']
                        bg[name] = f"{row['background_percent']:.1f}"
                        for key in ['common_retained_core_RNA', 'core_input_RNA', 'common_excluded_core_RNA', 'queries']:
                            if key in bg:
                                assert bg[key] == record[key], (region, key)
                            bg[key] = record[key]
                    else:
                        assert donor == 'LUNG' and index == 1 and name == '200_default', 'Unexpected failed or pending setting'
                        bg[name] = 'failed'
                    rows.append(row)
            background.append(bg)
    assert len(rows) == 32 and sum(r['returncode'] == 0 for r in rows) == 31
    output.mkdir(parents=True, exist_ok=True)
    (output / 'proseg_summary.json').write_text(json.dumps(rows, indent=2) + '\n')
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with (output / 'proseg_summary.csv').open('w', newline='') as target:
        writer = csv.DictWriter(target, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    content = [r'\begin{table}[ht]\centering\small',
        r'\begin{tabular}{lrrrrr}\toprule',
        r'Core&Retained/input RNA&\multicolumn{2}{c}{Default factorization}&\multicolumn{2}{c}{No factorization}\\',
        r'&&$200\,\mu\mathrm m$&$800\,\mu\mathrm m$&$200\,\mu\mathrm m$&$800\,\mu\mathrm m$\\\midrule']
    for row in background:
        content.append(f"{row['core']}&{row['common_retained_core_RNA']:,}/{row['core_input_RNA']:,}&{row['200_default']}&{row['800_default']}&{row['200_no_factorization']}&{row['800_no_factorization']}" + r'\\')
    content += [r'\bottomrule\end{tabular}',
        r'\caption{Percentage of retained core transcripts assigned to background. Each row uses the intersection of core transcript identifiers across completed settings; the second column reports its size against the original core input. Input width describes the surrounding window, while evaluation always uses the central core.}\label{tab:context-background}',
        r'\end{table}']
    content += [r'\begin{table}[ht]\centering\small',
        r'\begin{tabular}{llrrrr}\toprule',
        r'Core&Setting&Changed counts&Reach (\%)&With BG (\%)&Fixed RNA moves\\\midrule']
    for row in rows:
        core = f"{row['donor']} {int(row['region'].split('_')[-1])}"
        setting = f"{row['input_um']} / " + ('D' if row['factorization'] == 'default' else 'NF')
        if row['returncode'] != 0:
            content.append(f'{core}&{setting}' + r'&\multicolumn{4}{c}{Factorization failed}\\')
            continue
        content.append(f"{core}&{setting}&{row['changed_counts']:,}&{row['reach_percent']:.1f}&{row['background_related_percent']:.1f}&{row['changed_fixed_RNA']:,}" + r'\\')
    content += [r'\bottomrule\end{tabular}',
        r'\caption{All paired Proseg settings. D: default factorization; NF: no factorization. Reach is the fraction of changed cell--gene counts individually attainable within the primary atlas at $4\,\mu\mathrm m^2$. With BG counts changed entries touched by at least one transfer to or from background; these may also involve cell-to-cell transfers. Fixed RNA moves count reassigned transcripts outside editable ribbons. Every completed full assignment lies outside the family. Each core uses a fixed set of cell--gene questions across its completed settings.}\label{tab:context-reach}',
        r'\end{table}']
    content.append(r'\FloatBarrier')
    content.append(r'In the original lung cores, the matched no-factorization background fractions change from 42.5\% to 19.4\% and from 89.8\% to 62.1\% when surrounding tissue is included. Thus the extreme background rates in the small-crop comparison depend on the input context. Larger inputs still produce many changes outside the primary family (Table~\ref{tab:context-reach}). These results diagnose the comparison setting; they do not establish which assignments are anatomically correct or isolate context effects from stochastic variability.')
    content.append('')
    content.append(r'The public artifacts retain core transcript assignments for every completed run, the complete setting inventory, countwise comparisons and audit results. Replaying those extracts reproduces the tables; repeating Proseg itself is stochastic. Native transcript outputs were reconciled with the full output count matrices before the core extracts were saved.')
    tex.write_text('\n'.join(content) + '\n')
    print('All 32 settings recorded; 31 completed and audited, one factorization failure retained.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--analysis', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--tex', type=Path, required=True)
    args = parser.parse_args()
    build(args.analysis, args.output, args.tex)
