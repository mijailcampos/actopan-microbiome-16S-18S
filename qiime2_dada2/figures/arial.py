# Register Microsoft Arial (corefonts) with matplotlib and make it the only font used.
import glob, os, warnings
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
FONTDIR = os.environ.get('ARIAL_DIR', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'font', 'ttf'))
for f in glob.glob(os.path.join(FONTDIR, '*.TTF')):
    fm.fontManager.addfont(f)
plt.rcParams.update({'font.family': 'Arial', 'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial',
                     'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold', 'mathtext.sf': 'Arial',
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})
warnings.filterwarnings('error', message='.*(findfont|Glyph|missing from font).*')
