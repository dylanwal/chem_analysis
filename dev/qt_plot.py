"""
Demonstrate use of GLLinePlotItem to draw cross-sections of a surface.
"""
import pathlib

import numpy as np
import plotly.graph_objs as go
import pyarrow as pa


########################################################################
# pyarrow
def numpy_to_feather(arr: np.ndarray, file_path: str | pathlib.Path):
    if not isinstance(file_path, pathlib.Path):
        file_path = pathlib.Path(file_path)
    if file_path.suffix != ".feather":
        file_path = file_path.with_suffix(".feather")

    arrays = [pa.array(col) for col in arr]
    names = [str(i) for i in range(len(arrays))]
    batch = pa.RecordBatch.from_arrays(arrays, names=names)
    with pa.OSFile(str(file_path), 'wb') as sink:
        with pa.RecordBatchStreamWriter(sink, batch.schema) as writer:
            writer.write_batch(batch)


def feather_to_numpy(path):
    source = pa.memory_map(path, 'r')
    table = pa.ipc.RecordBatchStreamReader(source).read_all()
    data = np.empty((table.num_columns, table.num_rows))
    for col in range(table.num_columns):
        data[col, :] = table.column(str(col))
    return data


def pack_time_series(x: np.ndarray, time_: np.ndarray, z: np.array) -> np.ndarray:
    data = np.empty((len(time_) + 1, len(x) + 1), dtype=z.dtype)
    data[0, 0] = 0
    data[0, 1:] = x
    data[1:, 0] = time_
    data[1:, 1:] = z
    return data


def unpack_time_series(data: np.array) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = data[0, 1:]
    y = data[1:, 0]
    z = data[1:, 1:]
    return x, y, z


def merge_feather_files(paths):
    data = []
    for path in paths:
        data.append(feather_to_numpy(path))
    return np.concatenate(data)


########################################################################
# numpy functions
def get_slice(x: np.ndarray, start=None, end=None) \
        -> slice:
    if start is None:
        start_ = 0
    else:
        start_ = np.argmin(np.abs(x - start))
    if end is None:
        end_ = -1
    else:
        end_ = np.argmin(np.abs(x - end))

    return slice(start_, end_)


def normalize(data):
    return data / np.max(data)


def smooth_data(data, sigma: int = 10):
    from scipy.ndimage import gaussian_filter
    import copy

    smooth_data = np.empty_like(data.data)
    for row in range(data.data.shape[0]):
        smooth_data[row] = gaussian_filter(data.data[row], sigma)
    new = copy.deepcopy(data)
    data.data = smooth_data
    return data


def align_max(data):
    max_indices = np.argmax(data, axis=1)
    roll_amount = -max_indices + max_indices[0]
    for i in range(1, data.shape[0]):
        data[i] = np.roll(data[i], roll_amount[i])
    return data


########################################################################
# data
class IRData:
    def __init__(self, times, wavenumber, data):
        self.time = times
        self.time_zeroed = self.time - self.time[0]
        self.wavenumber = wavenumber
        self.data = data

    @property
    def x(self):
        return self.wavenumber

    @classmethod
    def from_file(cls, data):
        wavenumber = data[0, 1:]
        times = data[1:, 0]
        data = np.flip(data[1:, 1:], axis=1)

        data = data - np.mean(data[0:10, :], axis=0)
        return cls(times, wavenumber, data)


class NMRData:
    def __init__(self, times, ppm, data):
        self.time = times
        self.time_zeroed = self.time - self.time[0]
        self.ppm = ppm
        self.data = data

    @property
    def x(self):
        return self.ppm

    @classmethod
    def from_file(cls, data):
        ppm = data[0, 1:]
        times = data[1:, 0]
        data = data[1:, 1:]
        return cls(times, ppm, data)


########################################################################
# pymcr
from pymcr.constraints import ConstraintNonneg, Constraint
from pymcr.mcr import McrAR
from pymcr.regressors import OLS, NNLS
from pymcr.constraints import ConstraintNonneg, ConstraintNorm


class ConstraintConv(Constraint):
    """
    Conversion constraint. sum(C) = 1

    Parameters
    ----------
    copy : bool
        Make copy of input data, A; otherwise, overwrite (if mutable)
    """

    def __init__(self, index=None, copy=False):
        """ A must be non-negative"""
        super().__init__(copy)
        self.index = index

    def transform(self, A):
        """ Apply nonnegative constraint"""
        if self.index is not None:
            if self.copy:
                A = np.copy(A)
            selected_columns = A[:, self.index]
            row_sums = selected_columns.sum(axis=1)
            selected_columns /= row_sums[:, np.newaxis]
            A[:, self.index] = selected_columns
            return A

        row_sums = A.sum(axis=1)
        if self.copy:
            return A / row_sums[:, np.newaxis]
        else:
            A /= row_sums[:, np.newaxis]
            return A


class ConstraintNonneg(Constraint):
    """
    Non-negativity constraint. All negative entries made 0.

    Parameters
    ----------
    copy : bool
        Make copy of input data, A; otherwise, overwrite (if mutable)
    """

    def __init__(self, index=None, copy=False):
        """ A must be non-negative"""
        super().__init__(copy)
        self.index = index

    def transform(self, A):
        """ Apply nonnegative constraint"""
        if self.index is not None:
            if self.copy:
                A = np.copy(A)
            selected_columns = A[:, self.index]
            A[:, self.index] = selected_columns * (selected_columns > 0)
            return A

        if self.copy:
            return A * (A > 0)
        else:
            A *= (A > 0)
            return A


def plot_mca_results(mcrar, x, D, times, conv):
    # resulting spectra
    fig = go.Figure()
    for i in range(mcrar.ST_.shape[0]):
        fig.add_trace(go.Scatter(y=normalize(mcrar.ST_[i, :]), x=x, name=f"mca_{i}"))
    fig.add_trace(go.Scatter(y=normalize(D[0, :]), x=x, name="early"))
    fig.add_trace(go.Scatter(y=normalize(D[int(D.shape[0] / 2), :]), x=x, name="middle"))
    fig.add_trace(go.Scatter(y=normalize(D[-1, :]), x=x, name="late"))

    fig.update_layout(autosize=False, width=1200, height=600, font=dict(family="Arial", size=18, color="black"),
                      plot_bgcolor="white", showlegend=True)
    fig.update_xaxes(title="<b>wavenumber (cm-1) (min)</b>", tickprefix="<b>", ticksuffix="</b>", showline=True,
                     linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
                     gridwidth=1, gridcolor="lightgray", autorange="reversed")
    fig.update_yaxes(title="<b>normalized absorbance</b>", tickprefix="<b>", ticksuffix="</b>", showline=True,
                     linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
                     gridwidth=1, gridcolor="lightgray")
    fig.show()

    fig = go.Figure()
    for i in range(mcrar.C_.shape[1]):
        fig.add_trace(go.Scatter(x=times, y=mcrar.C_[:, i], name=f"mca_{i}"))
    fig.update_layout(autosize=False, width=800, height=600, font=dict(family="Arial", size=18, color="black"),
                      plot_bgcolor="white", showlegend=True)
    fig.update_xaxes(title="<b>rxn time (min)</b>", tickprefix="<b>", ticksuffix="</b>", showline=True,
                     linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
                     gridwidth=1, gridcolor="lightgray")
    fig.update_yaxes(title="<b>conversion</b>", tickprefix="<b>", ticksuffix="</b>", showline=True,
                     linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
                     gridwidth=1, gridcolor="lightgray", range=[0, 1])
    fig.show()

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=times, y=conv))
    fig.update_layout(autosize=False, width=800, height=600, font=dict(family="Arial", size=18, color="black"),
                      plot_bgcolor="white", showlegend=False)
    fig.update_xaxes(title="<b>rxn time (min)</b>", tickprefix="<b>", ticksuffix="</b>", showline=True,
                     linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
                     gridwidth=1, gridcolor="lightgray")
    fig.update_yaxes(title="<b>conversion</b>", tickprefix="<b>", ticksuffix="</b>", showline=True,
                     linewidth=5, mirror=True, linecolor='black', ticks="outside", tickwidth=4, showgrid=False,
                     gridwidth=1, gridcolor="lightgray", range=[0, 1])
    fig.show()



ir_data = IRData.from_file(
    merge_feather_files(
        [
            r"G:\Other computers\My Laptop\post_doc_2022\Data\polymerizations\DW2-7\DW2_7_flow_ATIR.feather",
            r"G:\Other computers\My Laptop\post_doc_2022\Data\polymerizations\DW2-7\DW2_7_flow_ATIR2.feather",
            r"G:\Other computers\My Laptop\post_doc_2022\Data\polymerizations\DW2-7\DW2_7_flow_ATIR3.feather"
        ]
    )
)


numpy_to_feather(np.column_stack((ir_data.x, ir_data.data[540, :])), "ir")


def plot():
    import numpy as np

    import pyqtgraph as pg
    from pyqtgraph.Qt import QtGui

    app = pg.mkQApp("Plotting Example")

    win = pg.GraphicsLayoutWidget(show=True, title="Basic plotting examples")
    win.resize(1000, 600)
    win.setWindowTitle('pyqtgraph example: Plotting')
    label = pg.LabelItem(justify='right')
    win.addItem(label)
    # Enable antialiasing for prettier plots
    pg.setConfigOptions(antialias=True)
    p1 = win.addPlot(title="Zoom on selected region", row=0, col=1)
    curve1 = p1.plot(ir_data.x[::10], ir_data.data[1, ::10], pen=(255,0,0), name="Red curve")

    p2 = win.addPlot(title="Time Selection", row=2, col=1)
    x= ir_data.time_zeroed[::10]
    y=np.ones_like(x)
    p2.plot(x, y, pen=(255,255,255,200))
    lr = pg.InfiniteLine(201, movable=True, bounds=(0, np.max(x)), label='x={value:0.2f}',
                         labelOpts={'position':0.1, 'color': (200,200,100), 'fill': (200,200,200,50), 'movable': True})
    p2.addItem(lr)
    p2.setLimits(xMin=np.min(y), xMax=np.max(y))


    def update_plot():
        num_curves = 3
        index = np.argmin(np.abs(ir_data.time_zeroed - lr.value()))
        p1.clearPlots()
        for i in range(1, num_curves + 1):
            color = pg.intColor(i, num_curves)
            curve = p1.plot(ir_data.x[::10], ir_data.data[i + index, ::10], pen=color, name=f"Curve {i}")
        p1.setLabels(left='Y-Axis', bottom='X-Axis')


    def update_pos():
        lr.setPos(lr.value())

    lr.sigPositionChanged.connect(update_plot)
    p2.sigXRangeChanged.connect(update_pos)
    update_plot()


    # cross hair
    vLine = pg.InfiniteLine(angle=90, movable=False)
    hLine = pg.InfiniteLine(angle=0, movable=False)
    p1.addItem(vLine, ignoreBounds=True)
    p1.addItem(hLine, ignoreBounds=True)
    vb = p1.vb

    def mouseMoved(evt):
        pos = evt
        if p1.sceneBoundingRect().contains(pos):
            mousePoint = vb.mapSceneToView(pos)
            index = int(mousePoint.x())
            label.setText(
                f"<span style='font-size: 12pt'>x={mousePoint.x()}</span>, "
                f"<span style='color: red'>y1={mousePoint.y()}</span>, " # data1[index]
                f"<span style='color: green'>y2={1}</span>") # data2[index]
            vLine.setPos(mousePoint.x())
            hLine.setPos(mousePoint.y())


    p1.scene().sigMouseMoved.connect(mouseMoved)

    spinBox = pg.SpinBox()
    win.addItem(spinBox)

    if __name__ == '__main__':
        pg.exec()


def plot2():
    """
    This example demonstrates the SpinBox widget, which is an extension of
    QDoubleSpinBox providing some advanced features:

      * SI-prefixed units
      * Non-linear stepping modes
      * Bounded/unbounded values

    """

    import pyqtgraph as pg
    from pyqtgraph.Qt import QtWidgets
    from pyqtgraph.parametertree import interact, ParameterTree, Parameter
    from pyqtgraph.dockarea.Dock import Dock
    from pyqtgraph.dockarea.DockArea import DockArea
    from pyqtgraph.Qt import QtCore

    app = pg.mkQApp("SpinBox Example")
    pg.setConfigOptions(antialias=True)
    win = QtWidgets.QMainWindow()
    area = DockArea()
    win.setCentralWidget(area)
    win.resize(1000, 500)
    win.setWindowTitle('pyqtgraph example: dockarea')

    d1 = Dock("Dock1", size=(800, 400), hideTitle=True)  ## give this dock the minimum possible size
    d2 = Dock("Dock2 - Console", size=(800, 100), hideTitle=True)
    d3 = Dock("Dock3", size=(200, 500), hideTitle=True)
    area.addDock(d1,'left')
    area.addDock(d2, 'bottom')  ## place d2 at right edge of dock area
    area.addDock(d3, 'right')  ## place d3 at bottom edge of d1


    class CustomViewBox(pg.ViewBox):
        def __init__(self, *args, **kwds):
            kwds['enableMenu'] = False
            pg.ViewBox.__init__(self, *args, **kwds)
            self.setMouseMode(self.RectMode)

        ## reimplement right-click to zoom out
        def mouseClickEvent(self, ev):
            if ev.button() == QtCore.Qt.MouseButton.RightButton:
                self.autoRange()

        ## reimplement mouseDragEvent to disable continuous axis zoom
        def mouseDragEvent(self, ev, axis=None):
            if axis is not None and ev.button() == QtCore.Qt.MouseButton.RightButton:
                ev.ignore()
            else:
                pg.ViewBox.mouseDragEvent(self, ev, axis=axis)

    vb = CustomViewBox()
    widget = pg.PlotWidget(viewBox=vb, enableMenu=False)
    d1.addWidget(widget)





    widget2 = pg.PlotWidget()
    x= ir_data.time_zeroed
    y=np.trapz(y=ir_data.data, axis=1)
    widget2.plot(x, y, pen=(255,255,255,200))
    lr = pg.InfiniteLine(201, movable=True, bounds=(0, np.max(x)), label='x={value:0.2f}',
                         labelOpts={'position':0.1, 'color': (200,200,100), 'fill': (200,200,200,50), 'movable': True})
    widget2.addItem(lr)
    widget2.setLimits(xMin=0, xMax=np.max(x))
    d2.addWidget(widget2)


    tree = ParameterTree()
    tree.setMinimumWidth(150)
    i3 = pg.TreeWidgetItem(["# Spectra"])
    i3_ = pg.SpinBox(value=5, int=True, dec=True, minStep=1, step=1, bounds=[1, 10])
    i3.setWidget(1, i3_)
    tree.addTopLevelItem(i3)
    d3.addWidget(tree)


    def update_plot():
        num_curves = i3_.value()
        index = np.argmin(np.abs(ir_data.time_zeroed - lr.value()))
        widget.clearPlots()
        for i in range(1, num_curves + 1):
            color = pg.intColor(i, num_curves)
            curve = widget.plot(ir_data.x[::10], ir_data.data[i + index, ::10], pen=color, name=f"Curve {i}")
        widget.setLabels(left='Y-Axis', bottom='X-Axis')


    lr.sigPositionChanged.connect(update_plot)
    i3_.sigValueChanged.connect(update_plot)
    update_plot()



    vLine = pg.InfiniteLine(movable=True)
    widget.addItem(vLine)
    i4 = pg.TreeWidgetItem(["x"])
    i4_ = pg.ValueLabel(formatStr='{value:4.0f}')
    i4.setWidget(1, i4_)
    tree.addTopLevelItem(i4)
    i5 = pg.TreeWidgetItem(["y"])
    i5_ = pg.ValueLabel(formatStr='{value:1.4f}')
    i5.setWidget(1, i5_)
    tree.addTopLevelItem(i5)

    def mouseMoved():
        x = vLine.value()
        i4_.setValue(x)
        index = np.argmin(np.abs(ir_data.x-x))
        y = ir_data.data[299, index]
        i5_.setValue(y)
        # label.setText(
        #     f"<span style='font-size: 12pt'>x={vLine.value()}</span>, "
        #     f"<span style='color: red'>y1={2}</span>, "  # data1[index]
        #     f"<span style='color: green'>y2={1}</span>")  # data2[index]
        #     # vLine.setPos(mousePoint.x())

    vLine.sigPositionChanged.connect(mouseMoved)


    win.show()
    if __name__ == '__main__':
        pg.exec()


if __name__ == '__main__':
    plot2()
