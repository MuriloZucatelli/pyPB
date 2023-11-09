import matplotlib as mpl

# mpl.use("agg")
import matplotlib.pyplot as plt
from numpy import sqrt
import seaborn as sns

# mpl.style.use('ggplot')


def set_plt_params(
    relative_fig_width=1.0, landscape=True, page_width=307.3, rescale_h=1
):
    fig_width_pt = page_width * relative_fig_width
    inches_per_pt = 1.0 / 72.27  # Convert pt to inch
    golden_mean = (sqrt(5.0) - 1.0) / 2.0  # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt  # width in inches

    if landscape:
        fig_height = fig_width * golden_mean  # height in inches
    else:
        fig_height = fig_width / golden_mean  # height in inches

    fig_height = fig_height * rescale_h
    fig_size = [fig_width, fig_height]
    params = {
        "font.family": "serif",
        "axes.labelsize": 7,
        "xtick.labelsize": 5,
        "ytick.labelsize": 5,
        "axes.labelcolor": "black",
        "ytick.color": "black",
        "xtick.color": "black",
        "legend.handlelength": 4,
        "legend.fontsize": 7,
        "legend.numpoints": 1,
        "lines.markersize": 3,
        # 'xtick.labelsize': 7,
        # 'ytick.labelsize': 7,
        "text.usetex": True,
        "figure.figsize": fig_size,
        "pgf.texsystem": "xelatex",
        "pgf.rcfonts": False,
    }

    plt.rcParams.update(params)


encoding = {
    "ff": ("fillet", "failed", "o"),
    "bf": ("box", "failed", "v"),
    "fn": ("fillet", "nonfailed", "o"),
    "bn": ("box", "nonfailed", "v"),
}


"""
    Para utilizar basta colocar no seu script:
    from plt_config import plt_config
    plt_config()
"""


def plt_config2(
    relative_fig_width=1.0, landscape=True, invert=False, page_width=210, rescale_h=1
):
    """_summary_

    Args:
        relative_fig_width (float, optional): _description_. Defaults to 0.8.
        landscape (bool, optional): _description_. Defaults to True.
        invert (bool, optional): _description_. Defaults to False.
        page_width (int, optional): _description_. Defaults to 210.
        rescale_h (int, optional): _description_. Defaults to 1.


        Avaliable styles for u flow pattern map:

        'Solarize_Light2', '_classic_test_patch', 'bmh', 'classic', 'dark_background', 'fast', 'fivethirtyeight',
        'ggplot','grayscale','seaborn','seaborn-bright','seaborn-colorblind', 'seaborn-dark', 'seaborn-dark-palette',
        'seaborn-darkgrid', 'seaborn-deep', 'seaborn-muted', 'seaborn-notebook', 'seaborn-paper', 'seaborn-pastel',
        'seaborn-poster','seaborn-talk','seaborn-ticks','seaborn-white','seaborn-whitegrid','tableau-colorblind10'

        Good styles for you case:
        'bmh', 'fast', 'ggplot', 'grayscale', 'seaborn-colorblind', 'seaborn-dark'

        see
        https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html
        for more information and aprreciation
    """

    # Define o renderizador:
    # mpl.use("QtAgg")
    # Agg rendering to a Tk canvas (requires TkInter). This backend can be activated in IPython with %matplotlib tk
    # Define o estilo primário:
    mpl.style.use("seaborn-colorblind")
    # mpl.style.use('fast')

    fig_width_pt = page_width * relative_fig_width
    mm_to_inches = 1.0 / 25.4  # Convert pt to mm
    golden_mean = (sqrt(5.0) - 1.0) / 2.0  # Aesthetic ratio
    fig_width = fig_width_pt * mm_to_inches  # width in mm

    if landscape:
        fig_height = fig_width * golden_mean  # height in mm
    else:
        fig_height = fig_width / golden_mean  # height in mm

    fig_height = fig_height * rescale_h

    fig_size = [fig_width, fig_height]

    if invert:
        fig_size = [fig_height*1.1, fig_width]

    # Modifica o estilo:
    params = {
        "font.family": "serif",
        "font.size": 12,
        "axes.labelsize": 11,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "grid.alpha": 0.3,
        "axes.labelcolor": "black",
        "ytick.color": "black",
        "xtick.color": "black",
        "legend.handlelength": 1.5,
        "legend.fontsize": 10,
        "legend.numpoints": 1,
        "legend.handletextpad": 0.2,
        "lines.markersize": 3,
        # 'xtick.labelsize': 7,
        # 'ytick.labelsize': 7,
        "text.usetex": True,
        "figure.figsize": fig_size,
        "pgf.texsystem": "xelatex",
        "pgf.rcfonts": False,
        "axes.formatter.min_exponent": 3,
    }

    plt.rcParams.update(params)

    encoding = {
        "ff": ("fillet", "failed", "o"),
        "bf": ("box", "failed", "v"),
        "fn": ("fillet", "nonfailed", "o"),
        "bn": ("box", "nonfailed", "v"),
    }
