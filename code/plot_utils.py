import string

def add_letters(ax, x=-0.08, y=1.02, fs=10, letter_offset=0):
    """
    adds bold letters a,b,c,... to the upper left corner of subplots
    :param ax: axis
    :param x: x location of text
    :param y: ylocation of text
    :param fs: fontsize
    :return:
    """
    letters = list(string.ascii_lowercase)
    try:
        ax.flat
        for il, tmp_ax in enumerate(ax.flat):
            tmp_ax.text(
                x,
                y,
                letters[il + letter_offset],
                weight="bold",
                horizontalalignment="center",
                verticalalignment="center",
                transform=tmp_ax.transAxes,
                fontsize=fs,
            )
    except AttributeError:
        ax.text(
            x,
            y,
            letters[letter_offset],
            weight="bold",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            fontsize=fs,
        )