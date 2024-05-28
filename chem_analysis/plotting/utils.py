from typing import Iterable


def rgb_string_to_tuple(rgb_string: str) -> tuple[int, ...]:
    """Convert an 'rgb(r,g,b)' string to a tuple of integers."""
    rgb_list = rgb_string.replace("rgb(", "").replace(")", "").split(",")
    return tuple(int(i) for i in rgb_list)


def hex_to_rgb(hex_color: str) -> tuple[int, ...]:
    """Convert hex color to RGB."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))


def rgb_to_hex(rgb_color: Iterable) -> str:
    """Convert RGB color to hex."""
    return '#{:02x}{:02x}{:02x}'.format(*rgb_color)


def darken_color(color: str, factor: float = 0.7) -> str:
    """
    Darken the given hex color by the specified factor.

    Args:
        color (str): The hex or rgb color to darken.
        factor (float): The factor by which to darken the color.
                        Should be between 0 and 1. Lower values result in darker colors.

    Returns:
        str: The darkened hex color.
    """
    if color.startswith("#"):
        rgb_color = hex_to_rgb(color)
    elif color.startswith("rgb"):
        rgb_color = rgb_string_to_tuple(color)
    else:
        raise ValueError("Not supported color format.")
    darkened_rgb = tuple(max(int(c * factor), 0) for c in rgb_color)
    return rgb_to_hex(darkened_rgb)
