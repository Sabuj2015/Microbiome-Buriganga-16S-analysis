#!/usr/bin/env python3
"""Assemble the multi-panel Figures 4, 5 and 6 from the regenerated panels.

Figure 4A (phylogenetic tree of dominant genera, drawn from the OTU representative sequences)
is included as an image (data/panels/).
"""

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

HERE = Path(__file__).resolve().parent
FIG = HERE.parent / "results" / "figures"
PANELS = HERE.parent / "data" / "panels"
Image.MAX_IMAGE_PIXELS = None


def font(size):
    for f in ("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", "DejaVuSans.ttf"):
        try:
            return ImageFont.truetype(f, size)
        except OSError:
            continue
    return ImageFont.load_default()


def label(im, text, xy=(20, 10)):
    ImageDraw.Draw(im).text(xy, text, fill="black", font=font(int(im.height * 0.035) + 20))
    return im


def open_rgb(path, height=None):
    im = Image.open(path).convert("RGB")
    if height:
        im = im.resize((round(im.width * height / im.height), height), Image.LANCZOS)
    return im


def side_by_side(images, gap=40):
    h = max(i.height for i in images)
    out = Image.new("RGB", (sum(i.width for i in images) + gap * (len(images) - 1), h), "white")
    x = 0
    for i in images:
        out.paste(i, (x, (h - i.height) // 2))
        x += i.width + gap
    return out


def stack(images, gap=40):
    w = max(i.width for i in images)
    out = Image.new("RGB", (w, sum(i.height for i in images) + gap * (len(images) - 1)), "white")
    y = 0
    for i in images:
        out.paste(i, ((w - i.width) // 2, y))
        y += i.height + gap
    return out


def save(im, name):
    im.save(FIG / f"{name}.tiff", dpi=(300, 300), compression="tiff_lzw")
    preview = im.copy()
    preview.thumbnail((2400, 2400))
    preview.save(FIG / f"{name}.png")


def main():
    b = label(open_rgb(FIG / "Figure_4B_reproduced.png", 1500), "B")
    c = label(open_rgb(FIG / "Figure_4C_reproduced.png", 1500), "C")
    a = open_rgb(PANELS / "Figure_4A_phylogenetic_tree.png", 3040)      # carries its own "A"
    save(side_by_side([a, stack([b, c])]), "Figure_4_reproduced")
    a5 = label(open_rgb(FIG / "Figure_5A_reproduced.png", 2700), "A")
    b5 = label(open_rgb(FIG / "Figure_5B_reproduced.png", 2700), "B")
    save(side_by_side([a5, b5]), "Figure_5_reproduced")
    a6 = label(open_rgb(FIG / "Figure_6A_reproduced.png", 2700), "A")
    b6 = label(open_rgb(FIG / "Figure_6B_reproduced.png", 2700), "B")
    save(side_by_side([a6, b6]), "Figure_6_reproduced")
    print("Figures 4-6 assembled")


if __name__ == "__main__":
    main()
