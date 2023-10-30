import pyqtgraph as pg

pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')


import sys
import numpy as np
import pyqtgraph as pg
from PyQt6.QtWidgets import QApplication, QMainWindow, QMenu, QVBoxLayout, QWidget, QStackedWidget
from PyQt6.QtCore import Qt
import PyQt6.QtGui as QtGui

import sys
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget

class MainMenu(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Chemistry Analysis Program")
        self.setGeometry(100, 100, 400, 300)

        central_widget = QWidget(self)
        self.setCentralWidget(central_widget)

        layout = QVBoxLayout()
        central_widget.setLayout(layout)

        self.nmr_button = QPushButton("NMR")
        # nmr_button.setFont(nmr_button.font().bold())
        self.nmr_button.setStyleSheet("font-size: 18px")
        self.nmr_button.clicked.connect(self.open_nmr)

        ir_button = QPushButton("IR")
        # ir_button.setFont(ir_button.font().bold())
        ir_button.setStyleSheet("font-size: 18px")
        ir_button.clicked.connect(self.open_ir)

        sec_button = QPushButton("SEC")
        # sec_button.setFont(sec_button.font().bold())
        sec_button.setStyleSheet("font-size: 18px")
        sec_button.clicked.connect(self.open_sec)

        layout.addWidget(self.nmr_button)
        layout.addWidget(ir_button)
        layout.addWidget(sec_button)

    def open_nmr(self):
        print("Opening NMR analysis")
        x = [0, 1, 2, 3, 4, 5]
        y = [0, 1, 0.5, 2, 1, 3]
        plot_widget = pg.PlotWidget()
        plot_widget.plot(x, y, pen=pg.mkPen("b"))

        layout = self.centralWidget().layout()
        layout.addWidget(plot_widget)
        self.nmr_button.hide()


    def open_ir(self):
        print("Opening IR analysis")

    def open_sec(self):
        print("Opening SEC analysis")

def main():
    app = QApplication(sys.argv)
    window = MainMenu()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()

# class MainMenu(QMainWindow):
#     def __init__(self):
#         super().__init__()
#
#         menubar = self.menuBar()
#         nmr_menu = menubar.addMenu("NMR")
#         ir_menu = menubar.addMenu("IR")
#         sec_menu = menubar.addMenu("SEC")
#
#         menu_style = "QMenu { font-size: 20px; }"
#         nmr_menu.setStyleSheet(menu_style)
#         ir_menu.setStyleSheet(menu_style)
#         sec_menu.setStyleSheet(menu_style)
#
#         nmr_action = QtGui.QAction("NMR Option", self)
#         nmr_action.triggered.connect(self.showNMRPlot)
#         ir_action = QtGui.QAction("IR Option", self)
#         sec_action = QtGui.QAction("SEC Option", self)
#
#         nmr_menu.addAction(nmr_action)
#         ir_menu.addAction(ir_action)
#         sec_menu.addAction(sec_action)
#
#         self.setWindowTitle("Main Menu Example")
#         self.setGeometry(100, 100, 600, 400)
#
#         self.stacked_widget = QStackedWidget(self)
#         self.setCentralWidget(self.stacked_widget)
#
#         self.nmr_plot_widget = self.createNMRPlot()
#         self.stacked_widget.addWidget(self.nmr_plot_widget)
#
#     def createNMRPlot(self):
#         plot_widget = pg.PlotWidget()
#         random_data = np.random.normal(size=100)
#         plot_widget.plot(random_data, pen=pg.mkPen('b'))
#         plot_widget.setTitle("Random NMR Data")
#         plot_widget.setLabel("left", "Amplitude")
#         plot_widget.setLabel("bottom", "Time")
#
#         central_widget = QWidget()
#         layout = QVBoxLayout()
#         layout.addWidget(plot_widget)
#         central_widget.setLayout(layout)
#
#         return central_widget
#
#     def showNMRPlot(self):
#         self.stacked_widget.setCurrentWidget(self.nmr_plot_widget)
#
#
# def main():
#     app = QApplication(sys.argv)
#     window = MainMenu()
#     window.show()
#     sys.exit(app.exec())
#
# if __name__ == "__main__":
#     main()
