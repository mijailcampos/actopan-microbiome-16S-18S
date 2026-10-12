"""Figure 1 — sampling sites, Actopan River mouth, Chachalacas, Veracruz, Mexico.

Inputs (same folder):
  sites.tsv                    site coordinates (WGS84), taken by the author in Google Earth
  data/mexico_states.geojson   Natural Earth 1:10m admin-1 (Mexico only)
  data/neighbors.geojson       Natural Earth 1:50m admin-0 (neighbouring countries)
Basemap: Esri World Imagery tiles (cached in data/tiles/).
Output: Figure1_sampling_map_600dpi.{png,tiff}
"""
import csv, io, json, math, os
import numpy as np
import requests
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import arial  # noqa: E402,F401
from matplotlib.patches import Polygon, Rectangle, FancyArrow
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'Figure1_sampling_map_600dpi')

# Map extent (WGS84)
LON_MIN, LON_MAX = -96.3318, -96.3192
LAT_MIN, LAT_MAX = 19.4108, 19.4215
ZOOM = 18
TILE_URL = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}'
ATTRIBUTION = 'Imagery: Esri, Maxar, Earthstar Geographics, and the GIS User Community'


# ---------- basemap ----------
def deg2num(lat, lon, z):
    n = 2 ** z
    x = (lon + 180) / 360 * n
    y = (1 - math.log(math.tan(math.radians(lat)) + 1 / math.cos(math.radians(lat))) / math.pi) / 2 * n
    return x, y


def num2deg(x, y, z):
    n = 2 ** z
    lon = x / n * 360 - 180
    lat = math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * y / n))))
    return lat, lon


def basemap():
    cache = os.path.join(HERE, 'data', 'tiles', str(ZOOM))
    os.makedirs(cache, exist_ok=True)
    x0, y0 = deg2num(LAT_MAX, LON_MIN, ZOOM)
    x1, y1 = deg2num(LAT_MIN, LON_MAX, ZOOM)
    tx0, ty0, tx1, ty1 = int(x0), int(y0), int(x1), int(y1)
    img = Image.new('RGB', (256 * (tx1 - tx0 + 1), 256 * (ty1 - ty0 + 1)))
    s = requests.Session()
    s.headers['User-Agent'] = 'scientific-figure/1.0'
    for tx in range(tx0, tx1 + 1):
        for ty in range(ty0, ty1 + 1):
            f = os.path.join(cache, f'{ty}_{tx}.jpg')
            if not os.path.exists(f):
                r = s.get(TILE_URL.format(z=ZOOM, y=ty, x=tx), timeout=30)
                r.raise_for_status()
                open(f, 'wb').write(r.content)
            img.paste(Image.open(f).convert('RGB'), (256 * (tx - tx0), 256 * (ty - ty0)))
    # crop to the exact extent (pixel space of the mosaic)
    crop = (round((x0 - tx0) * 256), round((y0 - ty0) * 256), round((x1 - tx0) * 256), round((y1 - ty0) * 256))
    img = img.crop(crop)
    # Web Mercator rows are not linear in latitude: resample rows to a linear-latitude grid
    arr = np.asarray(img)
    h = arr.shape[0]
    lats = np.linspace(LAT_MAX, LAT_MIN, h)
    rows = [(deg2num(la, LON_MIN, ZOOM)[1] - y0) * 256 for la in lats]
    rows = np.clip(np.round(rows).astype(int), 0, h - 1)
    return arr[rows]


# ---------- helpers ----------
def dms(v, pos, neg):
    hemi = pos if v >= 0 else neg
    v = abs(v)
    d = int(v)
    m = int((v - d) * 60)
    s = round(((v - d) * 60 - m) * 60)
    if s == 60:
        m, s = m + 1, 0
    return f'{d}°{m:02d}′{s:02d}″{hemi}'


def draw_geojson(ax, path, filt=None, **kw):
    for f in json.load(open(path))['features']:
        if filt and not filt(f):
            continue
        g = f['geometry']
        polys = g['coordinates'] if g['type'] == 'MultiPolygon' else [g['coordinates']]
        for poly in polys:
            ax.add_patch(Polygon(np.array(poly[0]), closed=True, **kw))


# ---------- main ----------
def main():
    sites = [r for r in csv.DictReader(open(os.path.join(HERE, 'sites.tsv')), delimiter='\t') if r['type'] != 'feces']
    kx = math.cos(math.radians((LAT_MIN + LAT_MAX) / 2))  # metres per degree lon / lat

    plt.rcParams.update({'font.size': 8})
    fig_w = 15 / 2.54  # 15 cm wide
    aspect = (LAT_MAX - LAT_MIN) / ((LON_MAX - LON_MIN) * kx)
    fig = plt.figure(figsize=(fig_w, fig_w * aspect * 0.86 + 0.55))
    ax = fig.add_axes([0.11, 0.09, 0.86, 0.89])

    ax.imshow(basemap(), extent=[LON_MIN, LON_MAX, LAT_MIN, LAT_MAX], interpolation='lanczos')
    ax.set_xlim(LON_MIN, LON_MAX)
    ax.set_ylim(LAT_MIN, LAT_MAX)
    ax.set_aspect(1 / kx)

    # graticule ticks every 10 arc-seconds
    step = 10 / 3600
    xt = np.arange(math.ceil(LON_MIN / step) * step, LON_MAX, step)
    yt = np.arange(math.ceil(LAT_MIN / step) * step, LAT_MAX, step)
    ax.set_xticks(xt)
    ax.set_xticklabels([dms(v, 'E', 'W') for v in xt])
    ax.set_yticks(yt)
    ax.set_yticklabels([dms(v, 'N', 'S') for v in yt], rotation=90, va='center')
    ax.tick_params(direction='out', length=3, labelsize=7)

    # sites
    style = {
        'water': dict(marker='o', s=70, c='#2b8cbe', label='River water'),
        'oyster': dict(marker='^', s=46, c='#fdae61', label=r'$\it{C.\ virginica}$ (oyster)'),
    }
    for r in sites:
        st = style[r['type']]
        ax.scatter(float(r['longitude']), float(r['latitude']), marker=st['marker'], s=st['s'],
                   c=st['c'], edgecolors='white', linewidths=1.0, zorder=5)
    # label offsets (points) chosen so labels do not overlap; OCW and O1–O10 share a location
    offsets = {'HZW': (-8, 7, 'right'), 'OCW': (-8, -11, 'right'), 'O1-O10': (8, -11, 'left'),
               'RMW': (8, 6, 'left'), 'GZW': (8, 6, 'left'), 'CZW': (-8, 6, 'right')}
    for r in sites:
        dx, dy, ha = offsets[r['site']]
        t = ax.annotate(r['label'], (float(r['longitude']), float(r['latitude'])), xytext=(dx, dy),
                        textcoords='offset points', ha=ha, va='center', fontsize=8, fontweight='bold',
                        color='white', zorder=6)
        t.set_path_effects([matplotlib.patheffects.withStroke(linewidth=2.2, foreground='black')])

    t = ax.text(-96.3245, 19.4199, 'Chachalacas\n(households)', ha='center', va='center', fontsize=8.5,
                fontweight='bold', color='white', zorder=6, linespacing=1.1)
    t.set_path_effects([matplotlib.patheffects.withStroke(linewidth=2.2, foreground='black')])

    leg = ax.legend(handles=[Line2D([], [], linestyle='none', markeredgecolor='white', markeredgewidth=1.0,
                                    marker=v['marker'], markersize=math.sqrt(v['s']), markerfacecolor=v['c'],
                                    label=v['label']) for v in style.values()],
                    loc='lower left', fontsize=7, frameon=True, framealpha=0.85, borderpad=0.6)
    leg.get_frame().set_edgecolor('none')

    # scale bar (200 m), bottom right
    m_per_deg_lon = 111320 * kx
    L = 200 / m_per_deg_lon
    x_end = LON_MAX - 0.07 * (LON_MAX - LON_MIN)
    y_sb = LAT_MIN + 0.06 * (LAT_MAX - LAT_MIN)
    hgt = 0.008 * (LAT_MAX - LAT_MIN)
    ax.add_patch(Rectangle((x_end - L - 0.02 * (LON_MAX - LON_MIN), y_sb - 3.2 * hgt), L + 0.05 * (LON_MAX - LON_MIN),
                           9 * hgt, facecolor='white', alpha=0.85, edgecolor='none', zorder=6))
    ax.add_patch(Rectangle((x_end - L, y_sb), L / 2, hgt, facecolor='black', zorder=7))
    ax.add_patch(Rectangle((x_end - L / 2, y_sb), L / 2, hgt, facecolor='white', edgecolor='black', lw=0.6, zorder=7))
    ax.add_patch(Rectangle((x_end - L, y_sb), L, hgt, facecolor='none', edgecolor='black', lw=0.6, zorder=7))
    for v, lab in [(x_end - L, '0'), (x_end - L / 2, '100'), (x_end, '200')]:
        ax.text(v, y_sb + 1.8 * hgt, lab, ha='center', va='bottom', fontsize=6.5, zorder=7)
    ax.text(x_end + 0.004 * (LON_MAX - LON_MIN), y_sb + hgt / 2, 'm', ha='left', va='center', fontsize=6.5, zorder=7)

    # north arrow, top right
    xa = LON_MAX - 0.05 * (LON_MAX - LON_MIN)
    ya = LAT_MAX - 0.17 * (LAT_MAX - LAT_MIN)
    ax.add_patch(FancyArrow(xa, ya, 0, 0.09 * (LAT_MAX - LAT_MIN), width=0.012 * (LON_MAX - LON_MIN),
                            head_width=0.035 * (LON_MAX - LON_MIN), head_length=0.04 * (LAT_MAX - LAT_MIN),
                            length_includes_head=True, color='white', ec='black', lw=0.6, zorder=7))
    t = ax.text(xa, ya - 0.035 * (LAT_MAX - LAT_MIN), 'N', ha='center', va='center', fontsize=9,
                fontweight='bold', color='white', zorder=7)
    t.set_path_effects([matplotlib.patheffects.withStroke(linewidth=2, foreground='black')])

    # attribution
    ax.text(0.995, 0.005, ATTRIBUTION, transform=ax.transAxes, ha='right', va='bottom', fontsize=4.8,
            color='white', zorder=7,
            path_effects=[matplotlib.patheffects.withStroke(linewidth=1.2, foreground='black')])

    # inset: Mexico with Veracruz highlighted
    ins = fig.add_axes([0.115, 0.70, 0.30, 0.275])
    ins.set_facecolor('#dfeef7')
    draw_geojson(ins, os.path.join(HERE, 'data', 'neighbors.geojson'), facecolor='#eeeeee', edgecolor='#999999', lw=0.3)
    draw_geojson(ins, os.path.join(HERE, 'data', 'mexico_states.geojson'),
                 facecolor='#fbfbfb', edgecolor='#7a7a7a', lw=0.25)
    draw_geojson(ins, os.path.join(HERE, 'data', 'mexico_states.geojson'),
                 filt=lambda f: f['properties'].get('name') == 'Veracruz', facecolor='#d7301f',
                 edgecolor='#7a7a7a', lw=0.25, alpha=0.75)
    lat_c = float(np.mean([float(r['latitude']) for r in sites]))
    lon_c = float(np.mean([float(r['longitude']) for r in sites]))
    ins.plot(lon_c, lat_c, marker='*', ms=7, mfc='yellow', mec='black', mew=0.5)
    ins.set_xlim(-118.5, -86.5)
    ins.set_ylim(14.3, 33.0)
    ins.set_aspect(1 / math.cos(math.radians(23)))
    ins.set_xticks([])
    ins.set_yticks([])
    for sp in ins.spines.values():
        sp.set_linewidth(0.6)
    ins.text(-94.3, 30.6, 'MEXICO', fontsize=6, ha='center', color='#444444')
    ins.text(-93.8, 22.0, 'Gulf of\nMexico', fontsize=5, ha='center', style='italic', color='#336699')

    for ext in ('png', 'tiff'):
        kw = {'pil_kwargs': {'compression': 'tiff_lzw'}} if ext == 'tiff' else {}
        fig.savefig(f'{OUT}.{ext}', dpi=600, bbox_inches='tight', pad_inches=0.03, **kw)
        print('Saved', f'{OUT}.{ext}')
    print(f'study-area centroid: {dms(lat_c, "N", "S")}, {dms(lon_c, "E", "W")} ({lat_c:.5f}, {lon_c:.5f})')


if __name__ == '__main__':
    import matplotlib.patheffects  # noqa: F401
    main()
