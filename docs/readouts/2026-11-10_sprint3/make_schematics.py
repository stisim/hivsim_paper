"""
Schematics for the sprint 3 readout (10 Nov 2026): how each model is decomposed into slots, and the swap plan.
Run from anywhere: python make_schematics.py. Writes PNG (200 dpi) and SVG to figures/.
Content from docs/swap_design.md, docs/status_2026-11-09.md and analyses 04-24 (as of 10 Nov).
"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

OUT = Path(__file__).resolve().parent / 'figures'
OUT.mkdir(exist_ok=True)

# Fixed colours for the readout (dataviz reference palette, slots 1-3; validated all-pairs there)
EMOD, OPTIMA, HIVSIM = '#2a78d6', '#eb6834', '#1baf7a'
INK, INK2, MUTED = '#0b0b0b', '#52514e', '#8a8984'
SURFACE, CELL, GRID = '#fcfcfb', '#f3f2ee', '#e1e0d9'
# Status palette (dataviz, fixed), always with an icon and a word
STATUS = {
    'done':        ('#0ca30c', '#ffffff', '✓', 'done'),
    'partly':      ('#fab219', '#0b0b0b', '◐', 'partly done'),
    'not started': ('#e1e0d9', '#0b0b0b', '○', 'not started'),
    'blocked':     ('#d03b3b', '#ffffff', '✕', 'blocked on a decision'),
}
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'svg.fonttype': 'none'})


def box(ax, x, y, w, h, fc, ec='none', lw=0, r=0.6, z=1):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f'round,pad=0,rounding_size={r}', fc=fc, ec=ec, lw=lw, zorder=z))


def canvas(w, h):
    fig = plt.figure(figsize=(w / 10, h / 10), facecolor=SURFACE)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, w)
    ax.set_ylim(0, h)
    ax.axis('off')
    return fig, ax


def save(fig, name):
    for ext in ('png', 'svg'):
        fig.savefig(OUT / f'{name}.{ext}', dpi=200, facecolor=SURFACE)
    plt.close(fig)


def status_chip(ax, x, y, key, size=9, word=True):
    """ Status icon in its colour, then the status word in ink """
    fc, tc, icon, label = STATUS[key]
    box(ax, x, y - 1.0, 2.0, 2.0, fc, ec=MUTED if key == 'not started' else 'none', lw=0.6, r=0.4, z=3)
    ax.text(x + 1.0, y, icon, ha='center', va='center', fontsize=size, color=tc, weight='bold', zorder=4)
    if word:
        ax.text(x + 2.6, y, label, ha='left', va='center', fontsize=size, color=INK, zorder=4)


# ---------------------------------------------------------------------------------------------------------------
# 1. Decomposition
# ---------------------------------------------------------------------------------------------------------------

def decomposition():
    W, H = 150, 100
    fig, ax = canvas(W, H)
    ax.text(1.5, 97, 'All three models are cut into the same slots; where code is fused, the slots move together',
            fontsize=15, weight='bold', color=INK, va='center')
    ax.text(1.5, 93.6, 'Coarse swap units today: C, V, RN, D.  Circled numbers are the shared interfaces (bottom).',
            fontsize=10.5, color=INK2, va='center')

    slots = [('D', 'Demography'), ('R', 'Risk, sex work'), ('N', 'Sexual contact'), ('T', 'Transmission'),
             ('H', 'Natural history\n+ ART biology'), ('S', 'Seeding'), ('C', 'Cascade +\nprevention'),
             ('M', 'Measurement'), ('X', 'Solver')]
    row = {k: i for i, (k, _) in enumerate(slots)}
    top, rh, gap = 71.5, 5.6, 0.5           # row i spans [top - (i+1)*rh + gap/2, top - i*rh - gap/2]
    def ry(i):
        return top - (i + 1) * rh + gap / 2
    x0, cw, cg = 17, 36, 1.2
    cols = [('EMOD′', EMOD), ('HIVsim', HIVSIM), ('Optima′', OPTIMA)]
    cx = [x0 + j * (cw + cg) for j in range(3)]

    # Slot labels
    for i, (k, name) in enumerate(slots):
        y = ry(i) + (rh - gap) / 2
        ax.text(1.5, y, k, fontsize=13, weight='bold', color=INK, va='center')
        ax.text(4.6, y, name, fontsize=9.5, color=INK2, va='center', linespacing=1.0)

    # Lineage strips and headers
    lineage = [
        ('EMOD binary →¹ EMOD′ →² EMOD′-fixed',
         ['¹ within run noise of EMOD (an. 04, 05)',
          '² Daniel’s 2 fixes: target score 108 → 7.5 (an. 17)']),
        ('HIVsim Malawi →¹ HIVsim on STIsim v4 (#613)',
         ['¹ slot factories = old make_sim: 0 arrays differ (an. 11)',
          '   + condom fix: prevalence 15–49, 2020: 8.1% (an. 16)']),
        ('Optima HIV 2.12.2 →¹ Optima′ ODE →² agents',
         ['¹ totals within 2.4% of Optima, 2000–20 (acceptance)',
          '² means within 0.35% of the ODE at 400k agents (an. 20)']),
    ]
    for j, ((name, color), (chain, notes)) in enumerate(zip(cols, lineage)):
        box(ax, cx[j], 79.2, cw, 11.6, CELL, r=0.6)
        ax.text(cx[j] + 1, 88.6, 'Lineage', fontsize=9, color=MUTED, va='center')
        ax.text(cx[j] + 1, 86.0, chain, fontsize=9.5, color=INK, weight='bold', va='center')
        for k, n in enumerate(notes):
            ax.text(cx[j] + 1, 83.3 - 2.3 * k, n, fontsize=9, color=INK2, va='center')
        box(ax, cx[j], 72.4, cw, 5.6, color, r=0.6)
        ax.text(cx[j] + cw / 2, 75.2, name, fontsize=13, weight='bold', color='#ffffff', ha='center', va='center')

    # Cells: (column, first slot, last slot, title, detail, interface tags)
    cells = [
        (0, 'D', 'D', 'EmodDemography', '1960.5 start, 1985 burn-in saved; 1%', ''),
        (0, 'R', 'R', 'EmodRisk', 'LOW/MED by node, HIGH = FSW or client', '③'),
        (0, 'N', 'N', 'EmodNetwork', 'pair queues, 4 relationship types', ''),
        (0, 'T', 'S', 'EmodHIV', 'per-act transmission; Weibull prognosis,\nCD4, ART mortality; seeding 1983\n'
                                 '(on STIsim people: adapters, an. 12)', ''),
        (0, 'C', 'C', 'Condensed cascade', "EMOD's table on EmodHost; fixed: re-enrol", '①②'),
        (0, 'M', 'M', 'Outputs (+ EMOD reports)', 'tracked source of each infection', '④'),
        (0, 'X', 'X', None, 'agents only', ''),
        (1, 'D', 'D', 'STIsim demographics', 'Malawi, from 1985, 10,000 agents', ''),
        (1, 'R', 'N', 'StructuredSexual', 'risk groups 0/1/2, FSW, client inside it;\n'
                                          'stable / casual / one-time / sex work', '③'),
        (1, 'T', 'S', 'sti.HIV  (+ sti.ART biology)', 'beta per act; acute, latent, falling CD4;\n'
                                                     'prevalence by group in 1985', ''),
        (1, 'C', 'C', 'HIVTest, ANCTest, sti.ART, VMMC', 'ART to a number; on StisimHost', '①'),
        (1, 'M', 'M', 'Outputs', 'tracked source of each infection', '④'),
        (1, 'X', 'X', None, 'agents only', ''),
        (2, 'D', 'R', '39 populations', 'births, deaths, ageing, risk transfers,\n'
                                        'ForcePopsize; FSW / Clients / MSM', '③'),
        (2, 'N', 'T', 'Mean-field force of infection', 'acts per year between populations;\nno partnerships, no edges', ''),
        (2, 'H', 'H', '6 CD4 bands', 'death by band; CD4 moves on ART', ''),
        (2, 'S', 'S', 'Initial state by population', '1990', ''),
        (2, 'C', 'C', "Optima's cascade table", 'rates, target, holds (an. 23); OptimaHost', '①②'),
        (2, 'M', 'M', 'Outputs(OptimaHost) / ODE flows', 'source drawn by stage share of the force', '④'),
        (2, 'X', 'X', 'ODE ↔ agents', 'optima_prime.py / agents.py', ''),
    ]
    for j, a, b, title, detail, tags in cells:
        y0 = ry(row[b])
        y1 = ry(row[a]) + rh - gap
        h = y1 - y0
        merged = a != b
        color = cols[j][1]
        box(ax, cx[j], y0, cw, h, '#ffffff' if title else SURFACE, ec=color if title else GRID,
            lw=1.6 if merged else 0.8, r=0.6, z=2)
        if merged:
            ax.text(cx[j] + cw - 1, y1 - 1.5, f'{a}–{b} fused' if row[b] - row[a] > 1 else f'{a}+{b} fused',
                    fontsize=9, color=INK2, ha='right', va='center', style='italic', zorder=3)
        mid = y0 + h / 2
        if title is None:
            ax.text(cx[j] + cw / 2, mid, detail, fontsize=9, color=MUTED, ha='center', va='center', zorder=3)
            continue
        nlines = detail.count('\n') + 1
        ty = mid + 0.95 * nlines if merged else mid + 1.15
        ax.text(cx[j] + 1.2, ty, title, fontsize=10.5, weight='bold', color=INK, va='center', zorder=3)
        ax.text(cx[j] + 1.2, ty - 1.4, detail, fontsize=9, color=INK2, va='top', zorder=3, linespacing=1.15)
        if tags:
            ax.text(cx[j] + cw - 1, ty, tags, fontsize=11, color=INK, ha='right', va='center', zorder=3)

    # Coarse swap units: brackets on the right
    bx = cx[2] + cw + 2.0
    groups = [('D', 'D', 'D', 'demography'), ('R', 'N', 'RN', 'risk + network'),
              ('T', 'S', 'V', 'HIV disease'), ('C', 'C', 'C', 'cascade')]
    for a, b, lab, name in groups:
        y0 = ry(row[b]) + 0.3
        y1 = ry(row[a]) + rh - gap - 0.3
        ax.plot([bx, bx + 1, bx + 1, bx], [y1, y1, y0, y0], color=INK2, lw=1.4, solid_capstyle='round')
        ax.text(bx + 2.2, (y0 + y1) / 2 + 0.9, lab, fontsize=13, weight='bold', color=INK, va='center')
        ax.text(bx + 2.2, (y0 + y1) / 2 - 1.1, name, fontsize=9, color=INK2, va='center')
    y0, y1 = ry(row['X']) + 0.3, ry(row['M']) + rh - gap - 0.3
    ax.plot([bx, bx + 1, bx + 1, bx], [y1, y1, y0, y0], color=MUTED, lw=1.0, ls=(0, (2, 1.5)))
    ax.text(bx + 2.2, (y0 + y1) / 2, 'own steps,\nOptima\npaths only\n\nOptima′: D+R,\nN+T cross\nunits', fontsize=9, color=INK2, va='center',
            linespacing=1.1)
    ax.text(bx, 75.2, 'Swap units', fontsize=10.5, weight='bold', color=INK, va='center')

    # Shared interfaces
    iy = 19.0
    ax.text(1.5, iy, 'Shared interfaces', fontsize=11.5, weight='bold', color=INK, va='center')
    ifaces = [
        ('① CascadeHost', 'EmodHost, StisimHost, OptimaHost (an. 22)', 'done'),
        ('② Cascade table', "EMOD's (an. 06), Optima's (an. 23); HIVsim's to write", 'partly'),
        ('③ Population label', "anyone's people → Optima's 39 (an. 22)", 'done'),
        ('④ Outputs', 'Fig 3A sources on every host (common/outputs.py)', 'done'),
        ('Contact-edge spec  (N → T)', "to write; Optima's T has no edges", 'blocked'),
        ('Infectiousness to Optima’s force', 'to write (H → T for agents)', 'blocked'),
    ]
    iw, ih = 47.3, 6.6
    for k, (name, where, st) in enumerate(ifaces):
        x = 1.5 + (k % 3) * (iw + 1.2)
        y = iy - 2.4 - ih - (k // 3) * (ih + 1.0)
        box(ax, x, y, iw, ih, CELL, r=0.6)
        ax.text(x + 1.2, y + ih - 1.9, name, fontsize=10.5, weight='bold', color=INK, va='center')
        ax.text(x + 1.2, y + 1.8, where, fontsize=9, color=INK2, va='center')
        fc, tc, icon, label = STATUS[st]
        label = {'done': 'built', 'partly': 'partly built', 'blocked': 'blocked'}[st]
        lx = x + iw - 1.2 - 2.0 - 0.6 * len(label) - 1.0
        box(ax, lx, y + ih - 2.9, 2.0, 2.0, fc, r=0.4, z=3)
        ax.text(lx + 1.0, y + ih - 1.9, icon, ha='center', va='center', fontsize=9, color=tc, weight='bold', zorder=4)
        ax.text(lx + 2.6, y + ih - 1.9, label, fontsize=9, color=INK, va='center')
    save(fig, 'decomposition')


# ---------------------------------------------------------------------------------------------------------------
# 2. Swap plan
# ---------------------------------------------------------------------------------------------------------------

def swap_plan():
    W, H = 150, 76
    fig, ax = canvas(W, H)
    ax.text(1.5, 73, 'Swap plan: the first path is run; on the core path only the cascade step is done',
            fontsize=15, weight='bold', color=INK, va='center')
    ax.text(1.5, 69.6, "Each box is one hybrid step (one coarse slot swapped), in the order C → V → RN → D. "
            "V and RN wait on how agents' partnerships meet Optima's mean-field transmission.",
            fontsize=10.5, color=INK2, va='center')

    names = {'C': 'cascade', 'V': 'HIV disease', 'RN': 'risk + network', 'D': 'demography',
             'MX': 'measurement +\nODE ↔ agents', 'V*': 'HIV disease', 'RN*': 'network + FOI', 'D*': 'demog. + risk'}
    paths = [
        ('HIVsim → EMOD′', 'first path, an. 16', ('HIVsim', HIVSIM), ('EMOD′-fixed', EMOD), '→',
         [('C', 'done'), ('V', 'done'), ('RN', 'done'), ('D', 'done')],
         ['All 5 hybrids run (1000 runs); both ends exact.',
          'Shares 2020: C +0.14 interrupted, +0.10 on ART;',
          'RN +0.21, V −0.11 on ART (cancel in part);',
          'D moves counts, not shares.']),
        ('EMOD′ ↔ Optima′', 'core path (Fig 3A)', ('EMOD′-fixed', EMOD), ('Optima′ ODE', OPTIMA), '↔',
         [('C', 'done'), ('V*', 'blocked'), ('RN*', 'blocked'), ('D*', 'blocked'), ('MX', 'partly')],
         ['C, both ways (an. 24): moves the shares 0.5–0.75 of the gap',
          '(interrupted: all of it), but not the total: EMOD’s cascade',
          'lowers infections in either host. The cascade carries the split.',
          'X and M checked on Optima′ alone (an. 20, 22).']),
        ('EMOD′ → HIVsim', 'same order, other end', ('EMOD′-fixed', EMOD), ('HIVsim', HIVSIM), '→',
         [('C', 'not started'), ('V', 'not started'), ('RN', 'not started'), ('D', 'not started')],
         ["Needs HIVsim's cascade as a table (D2) and HIVsim's",
          'V and RN on EMOD′’s people (adapters not written).']),
        ('HIVsim ↔ Optima′', '', ('HIVsim', HIVSIM), ('Optima′ ODE', OPTIMA), '↔',
         [('C', 'not started'), ('V*', 'blocked'), ('RN*', 'blocked'), ('D*', 'blocked'), ('MX', 'partly')],
         ["Optima's table runs on StisimHost already; not run.",
          'Same V/RN/D decision as the core path.']),
    ]
    rowh = 14.5
    ex, ew, sw, sg = 22.5, 11.5, 11.0, 1.0
    for i, (title, sub, start, end, arrow, steps, notes) in enumerate(paths):
        yc = 58.5 - i * rowh
        bh = 9.6
        if i == 1:
            box(ax, 0.8, yc - rowh / 2 + 0.6, W - 1.6, rowh - 1.2, '#f6f5f1', ec=GRID, lw=0.8, r=0.8, z=0)
        ax.text(1.8, yc + 1.4, title, fontsize=12, weight='bold', color=INK, va='center')
        if sub:
            ax.text(1.8, yc - 1.4, sub, fontsize=9.5, color=INK2, va='center')
        # endpoints
        x = ex
        box(ax, x, yc - bh / 2, ew, bh, start[1], r=0.6, z=2)
        ax.text(x + ew / 2, yc, start[0].replace('-fixed', '\n-fixed'), fontsize=10, weight='bold', color='#ffffff',
                ha='center', va='center', zorder=3, linespacing=1.0)
        x += ew + sg
        for key, st in steps:
            fc, tc, icon, label = STATUS[st]
            box(ax, x, yc - bh / 2, sw, bh, '#ffffff', ec=MUTED, lw=0.8, r=0.6, z=2)
            box(ax, x, yc + bh / 2 - 2.6, sw, 2.6, fc, ec=MUTED if st == 'not started' else 'none', lw=0.6, r=0.6, z=2.5)
            ax.text(x + sw / 2, yc + bh / 2 - 1.3, f'{icon} {label if st != "blocked" else "blocked"}',
                    fontsize=9, color=tc, weight='bold', ha='center', va='center', zorder=3)
            short = {'MX': 'M + X'}.get(key, key.rstrip('*'))
            ax.text(x + sw / 2, yc - 0.1, short, fontsize=13, weight='bold', color=INK, ha='center', va='center', zorder=3)
            ax.text(x + sw / 2, yc - 3.0, names[key], fontsize=9, color=INK2, ha='center', va='center', zorder=3,
                    linespacing=1.0)
            x += sw + sg
        box(ax, x, yc - bh / 2, ew, bh, end[1], r=0.6, z=2)
        ax.text(x + ew / 2, yc, end[0].replace('-fixed', '\n-fixed').replace(' ODE', '\nODE'), fontsize=10,
                weight='bold', color='#ffffff', ha='center', va='center', zorder=3, linespacing=1.0)
        ax.text((ex + x + ew) / 2, yc - bh / 2 - 1.2, arrow, fontsize=12, color=INK2, ha='center', va='center')
        # notes
        nx = ex + ew + sg + 5 * (sw + sg) + ew + 2.0
        for k, n in enumerate(notes):
            ax.text(nx, yc + 1.25 * (len(notes) - 1) - 2.5 * k, n, fontsize=9.5, color=INK, va='center')

    # Status legend
    ly = 4.0
    ax.text(1.8, ly, 'Status', fontsize=10, weight='bold', color=INK, va='center')
    lx = 12
    for key in ('done', 'partly', 'not started', 'blocked'):
        status_chip(ax, lx, ly, key, size=9.5)
        lx += 6 + 1.0 * len(STATUS[key][3]) + 4
    ax.text(W - 1.5, ly, 'V*, RN*, D*: Optima’s fused H, N+T, D+R', fontsize=9.5, color=INK2, ha='right', va='center')
    save(fig, 'swap_plan')


if __name__ == '__main__':
    decomposition()
    swap_plan()
    print('Wrote', *sorted(p.name for p in OUT.iterdir()))
