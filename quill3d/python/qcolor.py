"""
This module provides methods for the generation of QColor colormaps and color mixing.
"""

import numpy as np
import matplotlib.colors as mcol
from matplotlib import cm

__WEIGHT_RED, __WEIGHT_GREEN, __WEIGHT_BLUE = 0.2126, 0.7152, 0.0722

def linearize_colormap(cm, l=0.0, r=1.0, name='default'):
    """
    Converts colormap to a colormap with linear brightness profile.

    Parameters
    ----------
    cm
        colormap
    l
        brightness at the minimum value.
    r
        brightness at the maximum value.
    name
        the name of the new colormap
    Returns
    -------
    The converted colormap.

    """
    col_number = cm.N
    x = np.linspace(0.0, 1.0, col_number)
    expected_brightness = np.linspace(l, r, col_number)
    cm_values = cm(x).T[:3]
    new_cm_values = cm_values * expected_brightness / __pixel_brightness(cm_values)
    new_cm_values = np.apply_along_axis(__normalize_color, 0, new_cm_values)
    return mcol.LinearSegmentedColormap.from_list(name, colors=new_cm_values.T)


def create_qcolor_colormaps():
    """
    Creates QColor colormaps and registers them in matplotlib. The names are:
    qcolor_orange, qcolor_green, qcolor_blue, qcolor_purple, qcolor_red.
    """
    names = ['qcolor_orange', 'qcolor_green', 'qcolor_blue', 'qcolor_purple', 'qcolor_red']
    base_colors = ['#ff7f00', '#4daf4a', '#377eb8', '#984ea3', '#e41a1c'] # One of colorbrewer palettes

    for name, bc in zip(names, base_colors):
        tmp_map = mcol.LinearSegmentedColormap.from_list('tmp', ['white', bc], N=256)
        new_map = linearize_colormap(tmp_map, 1.0, 0.6, name)
        cm.register_cmap(name, new_map)


def mix_images(image_array):
    """
    Calculates a new image by mixing several images.

    Parameters
    ----------
    image_array
        Array of matplotlib.image.AxesImage

    Returns
    -------
    An array of the RGBA values of the resulting image.
    """
    def get_rgb_data(image):
        return image.cmap(image.norm(image.get_array()))

    rgb_array = map(get_rgb_data, image_array)

    # mixing of colors is based on x = (x1 + x2) / (1 + x1 x2) formula for 1 - R, 1 - G, 1 - B.
    from functools import reduce
    return reduce(lambda x, y: x * y / (2 - x - y + x * y), rgb_array)


def __pixel_brightness(rgb):
    """
    Calculates brightness based on RGB value according to sRGB luminance model.
    """
    value = __WEIGHT_RED * rgb[0] + __WEIGHT_GREEN * rgb[1] + __WEIGHT_BLUE * rgb[2]
    return value


def __normalize_color(rgb):
    """
    Normalizes rgb values so they don't exceed 1.0.
    """
    a = [__WEIGHT_RED, __WEIGHT_GREEN, __WEIGHT_BLUE]
    new_color = list(rgb)
    for i in range(3):
        if rgb[i] > 1.0:
            new_color[i] = 1.0
            j, k = {0, 1, 2} - {i}
            delta = a[i] * (rgb[i] - 1.0) / (a[j] + a[k])
            new_color[j] = rgb[j] + delta
            new_color[k] = rgb[k] + delta
            if new_color[j] > 1.0:
                new_color[k] += a[j] * (new_color[j] - 1) / a[k]
                new_color[j] = 1.0
                if new_color[k] > 1.0:
                    new_color[k] = 1.0
            elif new_color[k] > 1.0:
                new_color[j] += a[k] * (new_color[k] - 1) / a[j]
                new_color[k] = 1.0
                if new_color[j] > 1.0:
                    new_color[j] = 1.0
            return new_color
    return new_color


create_qcolor_colormaps()