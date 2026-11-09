"""Digitise Bansi-Matharu et al. 2025 Figure 3A (Malawi), EMOD and Optima HIV panels.

Each bar (one per year, 2000-2040) is sampled down its centre (3 px wide), every pixel
between the 100% and 0% ticks is assigned to the nearest of the four legend colours or
white, and the share of each category is its pixel count over the bar's coloured
pixels. 1 percentage point is about 1.9 px in this image, and boundary pixels are
blurred by the JPEG, so shares are good to about +-1-2 points; shares under ~1.5 are
anti-aliasing at the bar edges, not data. Checked against the values in the text: EMOD
2024 undiagnosed 33.3 (read 32.6), on ART 46.1 (45.5), 2040 undiagnosed 23.7 (23.0);
Optima 2024 on ART 34.8 (35.6), interrupted 8.4 (8.5).

Writes common/fig3a_malawi.csv. Run from the repo root: python analyses/18_fig3a_emod/digitise_fig3a.py
"""
import csv
import numpy as np
from PIL import Image

IMG = 'resources/bansi-matharu2025_fig3.jpg'
OUT = 'common/fig3a_malawi.csv'
# Pixel rows of the 100% and 0% ticks, found from the y-axis tick marks (x = 85-88)
PANELS = [('EMOD', 153, 340), ('Optima HIV', 405, 592)]
COLOURS = {'undiagnosed': (150, 210, 230), 'diagnosed_never_treated': (250, 190, 40),
           'on_art': (160, 200, 110), 'interrupted_art': (200, 90, 90), None: (255, 255, 255)}
CATS = [c for c in COLOURS if c]
HEADER = """\
# Bansi-Matharu et al. 2025 (Lancet Glob Health 13:e1006-19) Figure 3A, Malawi: share (%) of new
# infections by the source partner's treatment status, digitised by analyses/18_fig3a_emod/digitise_fig3a.py
# (+-1-2 points; shares under ~1.5 are edge anti-aliasing). Read with pd.read_csv(path, comment='#').
# Categories: undiagnosed; diagnosed_never_treated (diagnosed, ART-naive); on_art; interrupted_art
# (not receiving ART after having started).
# Caption totals, new infections in 2024 (national): EMOD 30582; Optima HIV 8459; Goals 10120; HIV Synthesis 22689.
"""


def bar_centres(im, row):
    coloured = im[row, 95:552].min(1) < 200
    centres, start = [], None
    for x, c in zip(range(95, 552), coloured):
        if c and start is None:
            start = x
        elif not c and start is not None:
            centres.append((start + x - 1) // 2)
            start = None
    return centres


def main():
    im = np.array(Image.open(IMG).convert('RGB')).astype(int)
    refs = np.array(list(COLOURS.values()))
    names = list(COLOURS)
    rows = []
    for model, y100, y0 in PANELS:
        xs = bar_centres(im, y0 - 5)
        assert len(xs) == 41, (model, len(xs))
        for year, x in zip(range(2000, 2041), xs):
            col = im[y100 - 3:y0 + 1, x - 1:x + 2].mean(1)
            lab = [names[i] for i in ((col[:, None, :] - refs[None]) ** 2).sum(2).argmin(1)]
            n = sum(l is not None for l in lab)
            rows += [(model, year, c, round(100 * lab.count(c) / n, 1)) for c in CATS]
    with open(OUT, 'w', newline='') as f:
        f.write(HEADER)
        w = csv.writer(f)
        w.writerow(['model', 'year', 'category', 'share'])
        w.writerows(rows)


if __name__ == '__main__':
    main()
