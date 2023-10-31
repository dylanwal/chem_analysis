import sys
import os

from PyQt6 import QtWidgets
from PyQt6 import QtCore
from PyQt6 import QtGui

from gui.qt.nmr import NMRView
from gui.qt.ir_array import IRArrayView


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, model=None):
        super().__init__()

        self.model = model

        self.width = 1000
        self.height = 600
        self.setWindowTitle("chem_analysis")
        self.settings = QtCore.QSettings('apps', 'settings')

        menu = MainMenu(self)
        self.setCentralWidget(menu)
        self.show()


class MainMenu(QtWidgets.QWidget):
    def __init__(self, window: QtWidgets.QMainWindow):
        super().__init__()
        self.window = window

        layout = QtWidgets.QVBoxLayout()
        self.setLayout(layout)

        nmr_button = QtWidgets.QPushButton("NMR")
        nmr_button.setStyleSheet("font-size: 18px")
        nmr_button.clicked.connect(self.open_nmr)

        ir_button = QtWidgets.QPushButton("IR")
        ir_button.setStyleSheet("font-size: 18px")
        ir_button.clicked.connect(self.open_ir)

        sec_button = QtWidgets.QPushButton("SEC")
        sec_button.setStyleSheet("font-size: 18px")
        sec_button.clicked.connect(self.open_sec)

        nmr_time_button = QtWidgets.QPushButton("NMR_time")
        nmr_time_button.setStyleSheet("font-size: 18px")
        nmr_time_button.clicked.connect(self.open_nmr_timeseries)

        ir_time_button = QtWidgets.QPushButton("IR_time")
        ir_time_button.setStyleSheet("font-size: 18px")
        ir_time_button.clicked.connect(self.open_ir_timeseries)

        sec_time_button = QtWidgets.QPushButton("SEC_time")
        sec_time_button.setStyleSheet("font-size: 18px")
        sec_time_button.clicked.connect(self.open_sec_timeseries)

        layout.addWidget(nmr_button)
        layout.addWidget(ir_button)
        layout.addWidget(sec_button)
        layout.addWidget(nmr_time_button)
        layout.addWidget(ir_time_button)
        layout.addWidget(sec_time_button)

    def open_nmr(self):
        print("Opening NMR analysis")
        self.window.setCentralWidget(NMRView(self.window))

    def open_ir(self):
        print("Opening IR analysis")

    def open_sec(self):
        print("Opening SEC analysis")

    def open_ir_timeseries(self):
        print("Opening NMR timeseries analysis")
        self.window.setCentralWidget(IRArrayView(self.window))

    def open_nmr_timeseries(self):
        print("Opening IR timeseries analysis")

    def open_sec_timeseries(self):
        print("Opening SEC timeseries analysis")


if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    window = MainWindow()
    app.exec()
