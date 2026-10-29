"""Plot the prespecified contiguous regions on the native breast cell mask.

Only the source raster and the frozen selection protocol are read. The figure
contains no transcript, interval, or completion information.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import zarr

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir', type=Path, required=True)
    parser.add_argument('--selection-protocol', type=Path,
                        default=HERE / 'selection_protocol.json')
    parser.add_argument('--source-provenance', type=Path,
                        default=HERE / 'source_provenance.json')
    parser.add_argument('--output-stem', type=Path,
                        default=HERE.parents[1] / 'manuscript/figures/contiguous_regions')
    args = parser.parse_args()
    protocol = json.loads(args.selection_protocol.read_text())
    provenance = json.loads(args.source_provenance.read_text())
    if protocol['source_archive_sha256'] != provenance['sha256']:
        raise ValueError('The source archive differs from the selection protocol')
    source = args.source_dir / 'cells.zarr.zip'
    expected = next(m['sha256'] for m in provenance['members']
                    if m['name'] == source.name)
    if hashlib.sha256(source.read_bytes()).hexdigest() != expected:
        raise ValueError('Cell-mask source checksum mismatch')
    with zarr.ZipStore(str(source), mode='r') as store:
        group = zarr.open_group(store=store, mode='r')
        mask = group['masks/1'][:]
        transform = group['masks/homogeneous_transform'][:]
    pixel = 1 / float(transform[0, 0])
    if list(mask.shape) != protocol['raster_shape_yx']:
        raise ValueError('Raster dimensions differ from the selection protocol')
    if pixel != protocol['pixel_size_um']:
        raise ValueError('Pixel scale differs from the selection protocol')

    # Thin grey interfaces preserve cell geometry without assigning colours to
    # particular cells. Downsample only the display, after finding interfaces.
    edge = np.zeros(mask.shape, dtype=bool)
    edge[:-1, :] |= mask[:-1, :] != mask[1:, :]
    edge[:, :-1] |= mask[:, :-1] != mask[:, 1:]
    shade = np.where(mask > 0, 0.935, 0.995).astype(np.float32)
    shade[edge] = 0.78
    display = shade[::2, ::2]
    height, width = mask.shape
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 11,
        'axes.titlesize': 12, 'axes.labelsize': 11.5,
        'axes.linewidth': 0.6, 'xtick.direction': 'out',
        'ytick.direction': 'out', 'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })
    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    ax.imshow(display, cmap='gray', vmin=0, vmax=1, origin='upper',
              extent=(0, width * pixel, height * pixel, 0),
              interpolation='nearest', rasterized=True)
    styles = {400: ('#0072B2', 'solid'), 600: ('#D55E00', (0, (6, 3)))}
    handles = []
    for region in sorted(protocol['regions'],
                         key=lambda r: r['nominal_side_um'], reverse=True):
        side = region['nominal_side_um']
        if side not in styles:
            raise ValueError(f'Unexpected prespecified region: {side}')
        x0, y0, x1, y1 = region['bbox_pixels']
        n = region['side_pixels']
        if [x0, y0, x1, y1] != [(width-n)//2, (height-n)//2,
                                (width-n)//2+n, (height-n)//2+n]:
            raise ValueError('A region is not the prespecified central square')
        color, linestyle = styles[side]
        ax.add_patch(Rectangle((x0*pixel, y0*pixel), (x1-x0)*pixel,
                               (y1-y0)*pixel, facecolor='none',
                               edgecolor=color, linestyle=linestyle,
                               linewidth=1.8, zorder=3))
        handles.append((side, Line2D([0], [0], color=color,
                                     linestyle=linestyle, linewidth=1.8,
                                     label=f'{side} µm region')))
    ax.legend(handles=[h for _, h in sorted(handles)], loc='lower left',
              bbox_to_anchor=(0, 1.015), ncol=2, frameon=False,
              borderaxespad=0, handlelength=3, columnspacing=2)
    ax.set_title('Prespecified contiguous regions', loc='left', pad=35)
    ax.set_xlabel('x (µm)')
    ax.set_ylabel('y (µm)')
    ax.set_xticks(np.arange(0, width*pixel, 200))
    ax.set_yticks(np.arange(0, height*pixel, 100))
    ax.set_xlim(0, width*pixel)
    ax.set_ylim(height*pixel, 0)
    ax.set_aspect('equal')
    ax.tick_params(length=3, width=0.6, colors='#333333')
    for spine in ax.spines.values():
        spine.set_color('#777777')
    fig.subplots_adjust(left=0.105, right=0.985, bottom=0.13, top=0.80)
    args.output_stem.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ('.pdf', '.png'):
        destination = args.output_stem.with_suffix(suffix)
        with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as temporary:
            temporary_path = Path(temporary.name)
        try:
            fig.savefig(temporary_path, dpi=300, bbox_inches='tight', pad_inches=0.08,
                        metadata={'Creator': 'Matplotlib'} if suffix == '.pdf' else None)
            os.replace(temporary_path, destination)
        finally:
            temporary_path.unlink(missing_ok=True)
        print(destination)
    plt.close(fig)


if __name__ == '__main__':
    main()
