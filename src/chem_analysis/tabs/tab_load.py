
from dash import html, dcc, dash_table



tab_load_layout = html.Div([
    html.H1('Load data'),
    dcc.Upload(
        id="upload-data",
        children=html.Div([
            'Drag and Drop or ',
            html.A('Select Files')
        ]),
        style={
            'width': '50%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '10px'
        },
        multiple=False,  # Allow multiple files to be uploaded
    ),
    html.Div(id='output-data-upload'),
    dash_table.DataTable(
        id='data_table',
        columns=[{"id": "time", "name": "time"}, {"id": "RI", "name": "RI"}],
        style_table={'overflowX': 'auto'},
        style_cell={
            'height': 'auto',
            # all three widths are needed
            'minWidth': '180px', 'width': '180px', 'maxWidth': '180px',
            'whiteSpace': 'normal'
        },
        persistence=True
                         ),
])
