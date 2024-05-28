from collections import OrderedDict

import plotly.graph_objs as go
import chem_analysis as ca

lib_path = r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library_color.json"
LIBRARY = ca.mass_spec.GCLibrary.from_JSON(lib_path)


COLORS = OrderedDict()
COLORS["misc"] = "545454"
COLORS["peroxide"] = "545400"
COLORS["alcohol TMS"] = "540000"
COLORS["alcohol acetylated"] = "633500"
COLORS["ketone"] = "046104"
COLORS["carboxylic acid TMS"] = "7FA8C7" #"000054"
COLORS["hydroxy acids TMS"] = "540054"
COLORS["hydroxy acids acetylated TMS"] = "540054"
COLORS["dicarboxylic acid TMS"] = "003030"


def lighter_shades(color: str, n: int) -> list[str]:
    hex_symbol_flag = False
    if color[0] == "#":
        hex_symbol_flag = True
        color = color[1:]

    # Convert hex color to RGB
    r = int(color[0:2], 16)
    g = int(color[2:4], 16)
    b = int(color[4:6], 16)

    # Calculate step size for each color channel
    r_step = (255 - r) / n
    g_step = (255 - g) / n
    b_step = (255 - b) / n

    # Generate lighter shades
    shades = []
    for i in range(n):
        r_new = int(r + i * r_step)
        g_new = int(g + i * g_step)
        b_new = int(b + i * b_step)
        new_color = "{:02x}{:02x}{:02x}".format(r_new, g_new, b_new)
        if hex_symbol_flag:
            new_color = "#" + new_color
        shades.append(new_color)

    return shades

def get_colors_from_scale(n: int = 10, color_scale: str = 'Jet') -> list[str]:
    import plotly.colors as colors
    color_scale = colors.get_colorscale(color_scale)
    return colors.sample_colorscale(color_scale, [i / (n - 1) for i in range(n)])


def color_square():
    fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
    n = 10
    count = 0
    for k, color in COLORS.items():
        fig.add_scatter(x=[count, count+1], y=[1, 1], name=k, mode="lines", fillcolor=f'#{color}', fill='tozeroy')
        count_2 = 0
        for color_ in lighter_shades(color, n):
            fig.add_scatter(x=[count, count+1], y=[1-(count_2)/n, 1-(count_2)/n],
                            showlegend=False, mode="lines", fillcolor=f'#{color_}', fill='tozeroy')
            count_2 += 1

        count += 1

    fig.show()


def main():
    with open(r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\colors.txt", 'r') as file:
        headers = file.readline()
        color_map = dict()
        for line in file.readlines():
            entries = line.split(",")
            color_map[entries[0]] = entries[1]

    for abbr, color in color_map.items():
        comp = LIBRARY.find_by_label(abbr)
        if comp is not None:
            comp.color = "#" + color


    # colors = get_colors_from_scale(len(LIBRARY.compounds))
    # for color, compound in zip(colors, LIBRARY.compounds):
    #     compound.color = color
    #
    LIBRARY.to_JSON(r"C:\Users\nicep\Desktop\research_wis\data\reference_data\gc_ms\decane\library_color.json",
                    binary=True, optimize=True, json_kwargs={"indent": 2}, overwrite=True)


if __name__ == "__main__":
    main()
