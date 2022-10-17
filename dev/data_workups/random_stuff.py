import plotly.graph_objects as go
from plotly.subplots import make_subplots

data_x1 = [0, 1, 2, 3, 4, 6, 10]
data_x2 = [0, 10, 20, 30, 40, 60, 100]
data_y1 = [0, 3555, 6114, 7395, 8016, 9130, 8646]
data_y2 = [1, 1.2, 1.18, 1.18, 1.22, 1.28, 1.52]

fig=make_subplots(specs=[[{"secondary_y": True}]])  # You need to create subplots figure with only one plot with the proper spec.
fig.update_layout(xaxis2= {'anchor': 'y', 'overlaying': 'x', 'side': 'top'})  # set the location and properties of the xaxis2

fig.add_trace(go.Scatter(x=data_x1, y=data_y1), secondary_y=False)
fig.add_trace(go.Scatter(x=data_x1, y=data_y2), secondary_y=True)
fig.add_trace(go.Scatter(x=data_x2, y=data_y1), secondary_y=False)
fig.add_trace(go.Scatter(x=data_x2, y=data_y2), secondary_y=True)

# for some reason reference to the x2 xaxis must be done after adding the trace
fig.data[2].update(xaxis='x2')
fig.data[3].update(xaxis='x2')
fig.show()