<<<<<<< HEAD
=======
<<<<<<< HEAD
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
    def __init__(self,mainWindow = None):
        super(self.__class__, self).__init__()
        self.mainWindow = mainWindow
        self.formLayout = QFormLayout()
        self.minDpLine = QLineEdit()
        self.minDpAsymLine = QLineEdit()
        self.maxDpAsymLine = QLineEdit()
        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)

        self.formLayout.addRow(self.tr("&MinDp"),self.minDpLine)
        self.formLayout.addRow(self.tr("&MinDpAsym"), self.minDpAsymLine)
        self.formLayout.addRow(self.tr("&MaxDpAsym"), self.maxDpAsymLine)

        self.buttonBox.accepted.connect(self.acceptInput)
        self.buttonBox.rejected.connect(self.reject)

        self.formLayout.addWidget(self.buttonBox)
        self.setLayout(self.formLayout)

    def acceptInput(self):
        try:
            self.minDp = float(self.minDpLine.text())
            self.minDpAsym = float(self.minDpAsymLine.text())
            self.maxDpAsym = float(self.maxDpAsymLine.text())
            self.accept()
        except:
            self.mainWindow.showError("Input data not valid.Please input again!")

    def getData(self):
        return self.minDp,self.minDpAsym, self.maxDpAsym


=======
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5
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
    def __init__(self,mainWindow = None):
        super(self.__class__, self).__init__()
        self.mainWindow = mainWindow
        self.formLayout = QFormLayout()
        self.minDpLine = QLineEdit()
        self.minDpAsymLine = QLineEdit()
        self.maxDpAsymLine = QLineEdit()
        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)

        self.formLayout.addRow(self.tr("&MinDp"),self.minDpLine)
        self.formLayout.addRow(self.tr("&MinDpAsym"), self.minDpAsymLine)
        self.formLayout.addRow(self.tr("&MaxDpAsym"), self.maxDpAsymLine)

        self.buttonBox.accepted.connect(self.acceptInput)
        self.buttonBox.rejected.connect(self.reject)

        self.formLayout.addWidget(self.buttonBox)
        self.setLayout(self.formLayout)

    def acceptInput(self):
        try:
            self.minDp = float(self.minDpLine.text())
            self.minDpAsym = float(self.minDpAsymLine.text())
            self.maxDpAsym = float(self.maxDpAsymLine.text())
            self.accept()
        except:
            self.mainWindow.show_error_dialog("Input processed_data not valid.Please input again!")

    def getData(self):
        return self.minDp,self.minDpAsym, self.maxDpAsym


<<<<<<< HEAD
=======
>>>>>>> ce5b242c238e196dfe20cf2038e62c260c5c95bd
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5
