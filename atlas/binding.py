"""Bind prepared-input certification to a previously checked physical atlas.

This is an integrity and applicability check, not a replacement for the rational
polygon validator. The saved geometry audit must name the exact content digest;
new/rebuilt atlases need analysis/verify_atlas.py before count certification.
An audit file is a reproducibility record, not a cryptographic proof against a
party that deliberately forges both inputs and validation records.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
from .transport import pixel_to_chart


def atlas_digest(atlas):
    body = {k: v for k, v in atlas.items() if k != 'atlas_id'}
    return hashlib.sha256(json.dumps(body, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def verify_prepared_atlas(atlas, source_dir, geometry_audit):
    """Reject altered/unbound atlases before any exact interval is calculated.

    The geometry audit is emitted by analysis/verify_atlas.py after exact mask,
    ribbon-overlap, membership and nuclear-fragment checks. Its atlas digest
    binds all patch definitions, observations, fixed counts and source hashes.
    """
    identity = atlas_digest(atlas)
    if atlas.get('atlas_id') != identity:
        raise ValueError('atlas content digest mismatch; rebuild and revalidate the atlas')
    if geometry_audit.get('failures') != 0:
        raise ValueError('a completed zero-failure geometry audit is required')
    records = [r for r in geometry_audit.get('coverage', [])
               if r.get('roi') == atlas.get('roi') and r.get('atlas_id') == identity]
    if len(records) != 1:
        raise ValueError('exact atlas identity absent or duplicated in geometry audit; run analysis/verify_atlas.py')
    source_dir = Path(source_dir).resolve()
    source_hashes = atlas.get('source_hashes', {})
    allowed = {'metadata.json', 'masks.npz', 'transcripts.parquet', 'transcripts.csv.gz', 'source_metadata.json'}
    if (not {'metadata.json', 'masks.npz'} <= set(source_hashes)
        or not set(source_hashes) <= allowed
        or len(set(source_hashes) & {'transcripts.parquet', 'transcripts.csv.gz'}) != 1):
        raise ValueError('metadata, mask and exactly one original transcript-file hash are required')
    for name, digest in source_hashes.items():
        path = source_dir / name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('atlas source-file hash mismatch: ' + name)
    pixel_size = F(atlas['pixel_size_um'])
    if pixel_size <= 0:
        raise ValueError('positive physical pixel size required')
    transform = json.loads((source_dir / 'metadata.json').read_text())['transform']
    if len(transform) != 4 or any(len(row) != 4 for row in transform):
        raise ValueError('expected explicit 4 by 4 source coordinate transform')
    tf = [[F(float(v)) for v in row] for row in transform]
    if (tf[0][0] <= 0 or tf[0][0] != tf[1][1]
        or any(tf[i][j] != 0 for i in (0, 1) for j in range(4) if j != i)
        or pixel_size != 1 / tf[0][0]):
        raise ValueError('this chart requires an isotropic unshifted XY source transform')
    cells = atlas['complete_cells']
    if len(cells) != len(set(cells)) or any(type(c) is not int or c <= 0 for c in cells):
        raise ValueError('unique positive complete-cell identifiers required')
    source = {}
    counts = defaultdict(Counter)
    for rec in atlas['records']:
        ident = rec['id']
        if not isinstance(ident, str) or not ident or ident in source:
            raise ValueError('unique nonempty retained transcript identifiers required')
        if type(rec['cell']) is not int or rec['cell'] < 0 or not isinstance(rec['gene'], str):
            raise ValueError('invalid retained transcript owner or gene')
        source[ident] = rec
        counts[str(rec['cell'])][rec['gene']] += 1
    expected = {cell: dict(genes) for cell, genes in counts.items()}
    supplied = {cell: {g: n for g, n in genes.items() if n != 0}
                for cell, genes in atlas['source_counts'].items() if any(genes.values())}
    if supplied != expected:
        raise ValueError('fixed source counts do not reconstruct from retained observations')
    patch_ids, moved_ids = set(), set()
    x0, y0 = atlas['bbox_pixels'][:2]
    for patch in atlas['patches']:
        pid = patch['id']
        if not isinstance(pid, str) or not pid or pid in patch_ids:
            raise ValueError('unique nonempty patch identifiers required')
        patch_ids.add(pid)
        labels = patch['cell_labels']
        if len(labels) != 2 or labels[0] == labels[1] or any(type(c) is not int or c < 0 for c in labels):
            raise ValueError('patch must separate two distinct valid labels')
        if any(c != 0 and c not in cells for c in labels):
            raise ValueError('patch changes an ineligible reference cell')
        if patch['sx'] not in (-1, 1) or patch['sy'] not in (-1, 1):
            raise ValueError('invalid physical chart orientation')
        width = F(patch['halfwidth_pixels'])
        if not 0 < F(patch['band_halfwidth_pixels']) < width:
            raise ValueError('strict fixed-rail displacement band required')
        if F(patch['area_scale_um2']) != pixel_size * pixel_size / 2:
            raise ValueError('physical Jacobian mismatch')
        lo, hi = map(F, patch['domain_pixels'])
        if lo != 0 or hi <= 0 or F(patch['slope']) <= 0 or F(patch['clearance_pixels']) <= 0:
            raise ValueError('invalid scalar patch domain, slope or clearance')
        for point in patch['points']:
            ident = point['transcript_id']
            if ident not in source or ident in moved_ids:
                raise ValueError('editable transcript is absent or belongs to multiple patches')
            moved_ids.add(ident)
            rec = source[ident]
            side = point['source']
            if side not in (0, 1) or point['gene'] != rec['gene'] or labels[side] != rec['cell']:
                raise ValueError('patch observation gene/source owner mismatch')
            # The ingestion contract retains the exact binary input coordinates.
            x = F(float(rec['x_um'])) / pixel_size - x0
            y = F(float(rec['y_um'])) / pixel_size - y0
            t, h = pixel_to_chart(patch, x, y)
            if (F(point['t']), F(point['x'])) != (t, h) or not lo < t < hi or not abs(h) < width:
                raise ValueError('patch observation coordinates/membership mismatch')
    return {'atlas_id': identity, 'roi': atlas['roi'],
            'patches': len(patch_ids), 'retained_transcripts': len(source),
            'editable_transcripts': len(moved_ids),
            'geometry_audit_binding': 'exact atlas content digest and original source hashes',
            'scope': 'one verified fixed-seam atlas; untouched boundary portions remain fixed'}
