import numpy as np
import plotly.graph_objects as go

def main():
    x = np.arange(11)
    y = np.array([0,1,2,3,4,5,4,3,2,1,0])
    y += -10

    print(np.trapezoid(x=x, y=y))
    print(np.trapezoid(x=x, y=y) - 0.5*(x[-1] - x[0]) *(y[0] + y[-1]))

    fig = go.Figure()
    fig.add_scatter(x=x, y=y, mode='markers')
    fig.show("browser")

if __name__ == "__main__":
    main()
