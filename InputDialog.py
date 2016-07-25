from InfoTableStyle import *
from CustomButton import *
import sys
import os
from PySide.QtGui import *
from PySide.QtCore import *
import settings
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import time


class InputForm(QDialog):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.formLayout = QFormLayout()
        self.minDp = QLineEdit()
        self.minDpAsym = QLineEdit()
        self.maxDpAsym = QLineEdit()
        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)

        self.formLayout.addRow(self.tr("&MinDp"),self.minDp)
        self.formLayout.addRow(self.tr("&MinDpAsym"), self.minDpAsym)
        self.formLayout.addRow(self.tr("&MaxDpAsym"), self.maxDpAsym)

        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)

        self.formLayout.addWidget(self.buttonBox)
        self.setLayout(self.formLayout)

    def getData(self):
        return self.minDp.text(),self.minDpAsym.text(), self.maxDpAsym.text()


