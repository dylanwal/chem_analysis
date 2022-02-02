
import pandas as pd
from dash.dependencies import Input, Output, State
from dash import html

from src.chem_analysis.server import app
from src.chem_analysis.analysis.utils.load_csv import load_csv


@app.callback([Output('output-data-upload', 'children'), Output('dataframe', 'data')],
              [Input('upload-data', 'contents')],
              [State('upload-data', 'filename')])
def upload_file(list_of_contents, filename):
    if list_of_contents is not None:
        result, df = load_csv(list_of_contents, filename)
        if df is not None:
            print(df.head())
            df = df.to_json(orient='split')
        return result, df
    return None, {}


@app.callback([Output('data_table', 'data'), Output('data_table', 'columns')],
              [Input('dataframe', 'data'), Input("data_table", "columns")])
def set_table(json_data, _):
    if json_data:
        print("hi")
        df = pd.read_json(json_data, orient='split')
        data = df.iloc[:10, :].to_dict('records')
        columns = [{"name": i, "id": i, "renamable": True, "deletable": True, "format": {"specifier": ".3f"}} for i
                   in df.columns]
        
        return data, columns
    return None, None




