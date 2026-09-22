#!/usr/bin/env python3
"""
Derive front/side overlap and detect crosshatch geometry from DJI image EXIF.

Usage:
    python3 overlap_from_exif.py                       # uses DEFAULT_FOLDER + DEFAULT_GSD below
    python3 overlap_from_exif.py "C:\\path\\to\\images" --gsd 0.75

  --gsd   ground sample distance in cm (from the Pix4D quality report)
  --sensor-w / --sensor-h  image dimensions in pixels (default 5280 x 3956)

Subfolders are searched recursively. Check the reported image count against the
"Total Images" figure in the matching quality report: if it is short, the flight
was split across sibling folders and you should point this at the parent.

Requires: pillow   (pip install pillow)
"""
import argparse, glob, math, os, re, sys
from PIL import Image

# ---------------------------------------------------------------- EXIF helpers

def _rational(v):
    try:
        return float(v)
    except Exception:
        return float(v[0]) / float(v[1])


def read_pose(path):
    """Return (lat, lon, yaw_deg or None) for one image."""
    lat = lon = yaw = None
    with open(path, "rb") as fh:
        head = fh.read(200_000)

    # XMP block: DJI writes FlightYawDegree, and often GPS as well.
    m = re.search(rb'drone-dji:FlightYawDegree="?([-+0-9.]+)', head)
    if m:
        yaw = float(m.group(1))
    m = re.search(rb'drone-dji:GpsLatitude="?([-+0-9.]+)', head)
    if m:
        lat = float(m.group(1))
    m = re.search(rb'drone-dji:GpsLongitude="?([-+0-9.]+)', head)
    if m:
        lon = float(m.group(1))

    if lat is None or lon is None:
        # Fall back to standard EXIF GPS IFD.
        img = Image.open(path)
        exif = img.getexif()
        gps = exif.get_ifd(0x8825)
        if gps:
            def dms(t, ref):
                d, mn, s = (_rational(x) for x in t)
                val = d + mn / 60 + s / 3600
                return -val if ref in ("S", "W") else val
            if 2 in gps and 4 in gps:
                lat = dms(gps[2], gps.get(1, "N"))
                lon = dms(gps[4], gps.get(3, "E"))
    return lat, lon, yaw


# ---------------------------------------------------------------- geometry

def to_metres(lats, lons):
    lat0 = sum(lats) / len(lats)
    lon0 = sum(lons) / len(lons)
    mlat = 111_132.9
    mlon = 111_320.0 * math.cos(math.radians(lat0))
    return [( (lo - lon0) * mlon, (la - lat0) * mlat ) for la, lo in zip(lats, lons)]


def dominant_headings(pts, yaws):
    """Heading of each exposure, from XMP yaw if present else from track."""
    if yaws and all(y is not None for y in yaws):
        return [y % 180.0 for y in yaws]
    h = []
    for i in range(len(pts)):
        j = min(i + 1, len(pts) - 1)
        k = max(i - 1, 0)
        dx = pts[j][0] - pts[k][0]
        dy = pts[j][1] - pts[k][1]
        h.append(math.degrees(math.atan2(dx, dy)) % 180.0)
    return h


def split_axes(headings, tol=30.0):
    """Group exposures into up to two perpendicular flight-line orientations."""
    hist = [0] * 180
    for h in headings:
        hist[int(h) % 180] += 1
    # circular smoothing over 180 deg
    sm = [sum(hist[(i + k) % 180] for k in range(-5, 6)) for i in range(180)]
    a1 = max(range(180), key=lambda i: sm[i])
    masked = [0 if min(abs(i - a1), 180 - abs(i - a1)) < 40 else sm[i] for i in range(180)]
    a2 = max(range(180), key=lambda i: masked[i]) if max(masked) > 0.15 * sm[a1] else None
    groups = {a1: [], a2: []} if a2 is not None else {a1: []}
    for idx, h in enumerate(headings):
        best = min(groups, key=lambda a: min(abs(h - a), 180 - abs(h - a)))
        if min(abs(h - best), 180 - abs(h - best)) <= tol:
            groups[best].append(idx)
    groups = {k: v for k, v in groups.items() if v}

    # Refine each axis to the circular mean (mod 180) of its members, so that
    # the projection below is aligned with the true flight-line direction.
    refined = {}
    for axis, idxs in groups.items():
        sx = sum(math.sin(math.radians(2 * headings[i])) for i in idxs)
        cy = sum(math.cos(math.radians(2 * headings[i])) for i in idxs)
        refined[(math.degrees(math.atan2(sx, cy)) / 2.0) % 180.0] = idxs
    return refined


def spacing_stats(pts, idxs, axis_deg):
    """Along-track and across-track spacing for one flight-line orientation."""
    th = math.radians(axis_deg)
    ux, uy = math.sin(th), math.cos(th)          # along flight line
    vx, vy = math.cos(th), -math.sin(th)         # across flight lines
    proj = [(pts[i][0] * ux + pts[i][1] * uy,
             pts[i][0] * vx + pts[i][1] * vy) for i in idxs]

    # cluster into lines by across-track coordinate
    proj.sort(key=lambda p: p[1])
    lines, cur = [], [proj[0]]
    for p in proj[1:]:
        if p[1] - cur[-1][1] > 8.0:              # 8 m gap starts a new line
            lines.append(cur); cur = [p]
        else:
            cur.append(p)
    lines.append(cur)
    lines = [l for l in lines if len(l) >= 3]

    along = []
    for l in lines:
        l.sort(key=lambda p: p[0])
        along += [l[i + 1][0] - l[i][0] for i in range(len(l) - 1)]
    along = [d for d in along if 0.1 < d < 100]
    centres = sorted(sum(p[1] for p in l) / len(l) for l in lines)
    across = [centres[i + 1] - centres[i] for i in range(len(centres) - 1)]

    med = lambda x: sorted(x)[len(x) // 2] if x else float("nan")
    return med(along), med(across), len(lines), len(idxs)


# ---------------------------------------------------------------- main

# DEFAULT_FOLDER = r"C:\Users\jakev\Downloads\TopSector1_23April2026+2026-04-23+\TopSector1_23April2026 2026-04-23\TopSector1_23April2026 2026-04-23 12_06_31 (UTC+02)\Remote-Control"
# DEFAULT_GSD = 0.88   # 23 April 2026, Top Sector, from the Pix4D quality report

DEFAULT_FOLDER = r"C:\Users\jakev\Downloads\23_March_2026_Top_Section_0.6\23_March_2026_Top_Section_0.6cm GSD\23_March_2026_Top_Section_0.6cm GSD 2026-03-23 15_14_13 (UTC+02)\Remote-Control"
DEFAULT_GSD = 0.75

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", nargs="?", default=DEFAULT_FOLDER)
    ap.add_argument("--gsd", type=float, default=DEFAULT_GSD, help="GSD in cm")
    ap.add_argument("--sensor-w", type=int, default=5280)
    ap.add_argument("--sensor-h", type=int, default=3956)
    a = ap.parse_args()

    exts = ("*.JPG", "*.jpg", "*.JPEG", "*.jpeg", "*.DNG", "*.dng")
    files = sorted(set(sum(
        [glob.glob(os.path.join(a.folder, "**", e), recursive=True) for e in exts], [])))
    if not files:
        sys.exit("No images found under " + a.folder)
    print(f"Folder: {a.folder}")
    print(f"Found {len(files)} image files (searched subfolders too)")

    lats, lons, yaws = [], [], []
    for f in files:
        la, lo, yw = read_pose(f)
        if la is None:
            continue
        lats.append(la); lons.append(lo); yaws.append(yw)
    print(f"{len(lats)} of {len(files)} images carried GPS")

    pts = to_metres(lats, lons)
    heads = dominant_headings(pts, yaws)
    groups = split_axes(heads)

    # footprint on the ground, metres
    fw = a.sensor_w * a.gsd / 100.0
    fh = a.sensor_h * a.gsd / 100.0
    print(f"GSD {a.gsd} cm -> footprint {fw:.1f} m x {fh:.1f} m\n")

    print(f"Flight-line orientations detected: {len(groups)}"
          f"  {'(CROSSHATCH)' if len(groups) > 1 else '(single grid)'}\n")

    for axis, idxs in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        al, ac, nl, ni = spacing_stats(pts, idxs, axis)
        # long sensor axis is normally across the flight line
        front_a = 1 - al / fh
        side_a  = 1 - ac / fw
        front_b = 1 - al / fw
        side_b  = 1 - ac / fh
        print(f"  axis {axis:5.1f} deg   {ni} exposures, {nl} lines")
        print(f"    along-track spacing  {al:6.2f} m")
        print(f"    across-track spacing {ac:6.2f} m")
        print(f"    if long axis across track : front {front_a*100:5.1f}%  side {side_a*100:5.1f}%"
              f"   product {(1-front_a)*(1-side_a):.4f}")
        print(f"    if long axis along track  : front {front_b*100:5.1f}%  side {side_b*100:5.1f}%"
              f"   product {(1-front_b)*(1-side_b):.4f}\n")


if __name__ == "__main__":
    main()