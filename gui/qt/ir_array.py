import pathlib

import numpy as np
from PyQt6 import QtWidgets
from PyQt6 import QtCore
from PyQt6 import QtGui

import pyqtgraph as pg
from pyqtgraph.parametertree import interact, ParameterTree, Parameter
from pyqtgraph.dockarea.Dock import Dock
from pyqtgraph.dockarea.DockArea import DockArea
from pyqtgraph.Qt import QtCore

from chem_analysis.utils.feather_format import unpack_and_merge_time_series_feather_files


pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')


class IRData:
    def __init__(self, wavenumber, times, data):
        self.time = times
        self.time_zeroed = self.time - self.time[0]
        self.wavenumber = wavenumber
        self.data = data

    @property
    def x(self):
        return self.wavenumber


def load_example() -> IRData:
    x, y, z = unpack_and_merge_time_series_feather_files(
            [
                r"G:\Other computers\My Laptop\post_doc_2022\Data\polymerizations\DW2-7\DW2_7_flow_ATIR.feather",
                r"G:\Other computers\My Laptop\post_doc_2022\Data\polymerizations\DW2-7\DW2_7_flow_ATIR2.feather",
                r"G:\Other computers\My Laptop\post_doc_2022\Data\polymerizations\DW2-7\DW2_7_flow_ATIR3.feather"
            ]
        )
    return IRData(x, y, z)


class IRArrayView(QtWidgets.QDockWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        IRMenuBar(self.parent().parent().window())
        self.stuff()

    def stuff(self):
        ir_data = load_example()

        area = DockArea()

        d1 = Dock("Dock1", size=(800, 400), hideTitle=True)  ## give this dock the minimum possible size
        d2 = Dock("Dock2 - Console", size=(800, 100), hideTitle=True)
        d3 = Dock("Dock3", size=(200, 500), hideTitle=True)
        area.addDock(d1, 'left')
        area.addDock(d2, 'bottom')  ## place d2 at right edge of dock area
        area.addDock(d3, 'right')  ## place d3 at bottom edge of d1
        self.layout().addWidget(area)

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
        x = ir_data.time_zeroed
        y = np.trapz(y=ir_data.data, axis=1)
        widget2.plot(x, y, pen=(255, 255, 255, 200))
        lr = pg.InfiniteLine(201, movable=True, bounds=(0, np.max(x)), label='x={value:0.2f}',
                             labelOpts={
                                 'position': 0.1, 'color': (200, 200, 100), 'fill': (200, 200, 200, 50), 'movable': True
                             })
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
            index = np.argmin(np.abs(ir_data.x - x))
            y = ir_data.data[299, index]
            i5_.setValue(y)
            # label.setText(
            #     f"<span style='font-size: 12pt'>x={vLine.value()}</span>, "
            #     f"<span style='color: red'>y1={2}</span>, "  # data1[index]
            #     f"<span style='color: green'>y2={1}</span>")  # data2[index]
            #     # vLine.setPos(mousePoint.x())

        vLine.sigPositionChanged.connect(mouseMoved)


class IRMenuBar:
    def __init__(self, window):
        menubar = window.menuBar()

        # main file
        menu_file = menubar.addMenu("File")
        action_open = menu_file.addAction("Open", self.open)
        action_save = menu_file.addAction("Save", self.save)
        action_quit = menu_file.addAction("Quit", self.quit)

        # processing
        menu_processing = menubar.addMenu("Processing")
        action_smooth = menu_processing.addAction("Smooth", self.smooth)
        action_baseline = menu_processing.addAction("Baseline Correction", self.baseline)

        # analysis
        menu_analysis = menubar.addMenu("Analysis")
        action_integrate = menu_analysis.addAction("Integrate", self.integrate)

        # Toolbars
        TBfile = window.addToolBar("File")
        TBfile.addAction(action_open)
        TBfile.addAction(action_save)

        open_icon = window.style().standardIcon(QtWidgets.QStyle.StandardPixmap.SP_DirOpenIcon)
        save_icon = window.style().standardIcon(QtWidgets.QStyle.StandardPixmap.SP_DialogSaveButton)

        action_open.setIcon(open_icon)
        action_save.setIcon(save_icon)

        TBfile.setAllowedAreas(QtCore.Qt.ToolBarArea.AllToolBarAreas)

        TBviewerNavigation = window.addToolBar("Viewer Navigation")

        # action_previous = QtGui.QAction("<", self)
        # action_next = QtGui.QAction(">", self)

        # TBviewerNavigation.addAction(action_previous)
        # TBviewerNavigation.addAction(action_next)
        TBviewerNavigation.setAllowedAreas(
            QtCore.Qt.ToolBarArea.TopToolBarArea | QtCore.Qt.ToolBarArea.BottomToolBarArea
        )

    def open(self):
        print("open")

    def save(self):
        print("save")

    def quit(self):
        print("quit")

    def smooth(self):
        print("smooth")

    def baseline(self):
        print("baseline")

    def integrate(self):
        print("integrate")
