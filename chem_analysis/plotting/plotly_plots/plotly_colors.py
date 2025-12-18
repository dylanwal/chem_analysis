
import plotly.colors as colors

def get_n_colors(num_colors: int = 1, style: str | None = None) -> list[str]:
    """Gets color for 2D plots."""
    color_scale = colors.get_colorscale(style or "phase")
    return colors.sample_colorscale(color_scale, [i / (num_colors - 1) for i in range(num_colors)])
