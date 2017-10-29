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


class KappaVarConfirmDialog(QDialog):
    def __init__(self,sigma,temp,dd1,iKappa1,dd2,iKapp2,solu,mainWindow = None):
        super(self.__class__, self).__init__()
        self.mainWindow = mainWindow
        self.formLayout = QFormLayout()
        self.sigmaLine = QLineEdit(str(sigma))
        self.tempLine = QLineEdit(str(temp))
        self.dd1Line = QLineEdit(str(dd1))
        self.iKappa1Line = QLineEdit(str(iKappa1))
        self.dd2Line = QLineEdit(str(dd2))
        self.iKappa2Line = QLineEdit(str(iKapp2))
        self.soluLine = QLineEdit(str(solu))
        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.formLayout.addRow(self.tr("&Sigma"),self.sigmaLine)
        self.formLayout.addRow(self.tr("&Temperature"), self.tempLine)
        self.formLayout.addRow(self.tr("&dry diameter(1)"), self.dd1Line)
        self.formLayout.addRow(self.tr("&iKappa(1)"), self.iKappa1Line)
        self.formLayout.addRow(self.tr("&dry diameter(2)"), self.dd2Line)
        self.formLayout.addRow(self.tr("&iKappa(2)"), self.iKappa2Line)
        self.formLayout.addRow(self.tr("&solubility"), self.soluLine)

        self.buttonBox.accepted.connect(self.acceptInput)
        self.buttonBox.rejected.connect(self.reject)

        self.formLayout.addWidget(self.buttonBox)
        self.setLayout(self.formLayout)

    def acceptInput(self):
        try:
            self.sigma = float(self.sigmaLine.text())
            self.temp = float(self.tempLine.text())
            self.dd1 = float(self.dd1Line.text())
            self.iKappa1 = float(self.iKappa1Line.text())
            self.dd2 = float(self.dd2Line.text())
            self.iKappa2 = float(self.iKappa2Line.text())
            self.solu = float(self.soluLine.text())
            self.accept()
        except:
            self.mainWindow.showError("Input data not valid.Please input again!")

    def getData(self):
        return self.sigma, self.temp, self.dd1, self.dd2, self.iKappa1, self.iKappa2, self.solu


=======
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


class KappaVarConfirmDialog(QDialog):
    def __init__(self,sigma,temp,dd1,iKappa1,dd2,iKapp2,solu,mainWindow = None):
        super(self.__class__, self).__init__()
        self.mainWindow = mainWindow
        self.formLayout = QFormLayout()
        self.sigmaLine = QLineEdit(str(sigma))
        self.tempLine = QLineEdit(str(temp))
        self.dd1Line = QLineEdit(str(dd1))
        self.iKappa1Line = QLineEdit(str(iKappa1))
        self.dd2Line = QLineEdit(str(dd2))
        self.iKappa2Line = QLineEdit(str(iKapp2))
        self.soluLine = QLineEdit(str(solu))
        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.formLayout.addRow(self.tr("&Sigma"),self.sigmaLine)
        self.formLayout.addRow(self.tr("&Temperature"), self.tempLine)
        self.formLayout.addRow(self.tr("&dry diameter(1)"), self.dd1Line)
        self.formLayout.addRow(self.tr("&iKappa(1)"), self.iKappa1Line)
        self.formLayout.addRow(self.tr("&dry diameter(2)"), self.dd2Line)
        self.formLayout.addRow(self.tr("&iKappa(2)"), self.iKappa2Line)
        self.formLayout.addRow(self.tr("&solubility"), self.soluLine)

        self.buttonBox.accepted.connect(self.acceptInput)
        self.buttonBox.rejected.connect(self.reject)

        self.formLayout.addWidget(self.buttonBox)
        self.setLayout(self.formLayout)

    def acceptInput(self):
        try:
            self.sigma = float(self.sigmaLine.text())
            self.temp = float(self.tempLine.text())
            self.dd1 = float(self.dd1Line.text())
            self.iKappa1 = float(self.iKappa1Line.text())
            self.dd2 = float(self.dd2Line.text())
            self.iKappa2 = float(self.iKappa2Line.text())
            self.solu = float(self.soluLine.text())
            self.accept()
        except:
            self.mainWindow.show_error_dialog("Input processed_data not valid.Please input again!")

    def getData(self):
        return self.sigma, self.temp, self.dd1, self.dd2, self.iKappa1, self.iKappa2, self.solu


>>>>>>> 8a1a2483f093860c8fd1098754e4b83032fc1d20
