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
from InputDialog import *
from KappaVarConfirmDialog import *

class PeakTextDataWidget(QWidget):
    def __init__(self, mainWindow=None):
        super(self.__class__,self).__init__(mainWindow)
        self.layout = QVBoxLayout()  #Vertical layout
        self.setLayout(self.layout)
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.infoTable = PeakDataTable(mainWindow)
        self.layout.addWidget(self.infoTable)

    def resize(self, parentWidth, parentHeight):
        self.setFixedWidth(parentWidth / 4)
        self.setFixedHeight(parentHeight)
        self.infoTable.resize(self.width(), self.height())

class PeakDataTable(QTableWidget):
    def resize(self, parentWidth, parentHeight):
        self.setFixedHeight(parentHeight)
        self.setFixedWidth(parentWidth)
        self.verticalHeader().setDefaultSectionSize(self.height() / 25)

    def __init__(self, mainWindow = None):
        QTableWidget.__init__(self)
        self.setColumnCount(1)
        self.setRowCount(0)
        self.setShowGrid(False)
        self.verticalHeader().setVisible(False)
        self.horizontalHeader().setVisible(False)
        self.horizontalHeader().setStretchLastSection(True)
        self.verticalHeader().setDefaultSectionSize(self.height() / 25)
        self.setWordWrap(True)
        self.resizeRowsToContents()
        self.setFrameStyle(QFrame.NoFrame)
        self.mainWindow = mainWindow

        #set background color
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, settings.infoAreaBackgroundColor)
        self.setPalette(palette)

    def updateGeneralInfo(self):
        header = TableHeader("Data Information")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        self.addMessage("Date",self.mainWindow.controller.date)
        self.addMessage("Time frame",self.mainWindow.controller.startTimeEntries[0] + " to " + self.mainWindow.controller.endTimeEntries[-1])
        self.addMessage("Time per run",self.mainWindow.controller.timeFrame)
        self.addMessage("Total run",self.mainWindow.controller.maxPeak)
        self.addMessage("CPC", self.mainWindow.controller.flowRate)

    def updateBasicPeakInfo(self):
        if self.rowCount() > 6:
            for i in range(self.rowCount() - 6):
                self.removeRow(self.rowCount()-1)
        header = TableHeader("Basic Peak Information")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        currPeak = self.mainWindow.controller.currPeak
        if self.mainWindow.controller.minPosCCNCList[currPeak] and self.mainWindow.controller.minPosCCNCList[currPeak]:
            self.addMessage("Status", "Valid for curve fit")
        else:
            self.addMessage("Status", "Invalid for curve fit", color='#EF5350')
        self.addMessage("Current run",currPeak + 1)
        self.addMessage("Saturation",self.mainWindow.controller.superSaturation)

    def updateSigFitPeakInfo(self):
        if self.rowCount() > 6:
            for i in range(self.rowCount() - 6):
                self.removeRow(self.rowCount() - 1)

        # Add basic information data
        header = TableHeader("Basic Peak Information")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        currPeak = self.mainWindow.controller.currPeak
        if self.mainWindow.controller.minPosCCNCList[currPeak] and self.mainWindow.controller.minPosCCNCList[currPeak]:
            self.addMessage("Status", "Valid for curve fit")
        else:
            self.addMessage("Status", "Invalid for curve fit",  color='#EF5350')
        self.addMessage("Current run", currPeak + 1)
        self.addMessage("Saturation", self.mainWindow.controller.superSaturation)

        # Add advance information data
        header = TableHeader("Advance Peak Information")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        if self.mainWindow.controller.usableForKappaCalList[currPeak]:
            self.addMessage("Status", "Valid for Kappa Cal")
        else:
            self.addMessage("Status", "Invalid for Kappa Cal",  color='#EF5350')
        self.addMessage('minDp',self.mainWindow.controller.minDp)
        self.addMessage('minDpAsym', self.mainWindow.controller.minDpAsym)
        self.addMessage('maxDpAsym', self.mainWindow.controller.maxDpAsym)
        self.addMessage("dp50", self.mainWindow.controller.dp50)
        self.addMessage("<dp50 counts",self.mainWindow.controller.dp50LessCount)
        self.addMessage(">dp50 counts",self.mainWindow.controller.dp50MoreCount)
        self.addMessage("dp50(Wet)",self.mainWindow.controller.dp50Wet)
        self.addMessage("dp50+20",self.mainWindow.controller.dp50Plus20)


    def addMessage(self, field,message, color = None):
        ##### add code here to process the message before printing
        message = str(message)
        item = TableItem(field, message, color)
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount()-1,0,item)

class PeakGraphWidget(QWidget):
    def resize(self, parentWidth, parentHeight):
        self.setFixedWidth(parentWidth * 3 / 4)
        self.setFixedHeight(parentHeight)
        self.tempAndMinView.resize(self.width(), self.height())
        self.controlWidget.resize(self.width(), self.height())
        self.dpAndDnlogView.resize(self.width(), self.height())

    def __init__(self, mainWindow=None):
        super(self.__class__, self).__init__(mainWindow)
        self.mainWindow = mainWindow
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.tempAndMinView = PeakRectFigureCanvas()
        self.dpAndDnlogView = PeakDpDnLogGraphWidget()
        self.controlWidget = PeakControlTabWidget(mainWindow)
        self.layout.addWidget(self.dpAndDnlogView)
        self.layout.addWidget(self.controlWidget)
        self.layout.addWidget(self.tempAndMinView)

class PeakControlTabWidget(QWidget):
    def resize(self, parentWidth, parentHeight):
        self.setFixedHeight(parentHeight * 1 / 10)
        self.setFixedWidth(parentWidth)
        self.nextButton.resize(self.width(),self.height())
        self.previousButton.resize(self.width(),self.height())
        self.removePeak.resize(self.width(), self.height())
        self.optimizeButton.resize(self.width(), self.height())
        self.updateSigData.resize(self.width(), self.height())
        self.calKappaButton.resize(self.width(), self.height())
        self.addSecond.resize(self.width(), self.height())
        self.subSecond.resize(self.width(), self.height())

    def __init__(self, mainWindow = None):
        self.mainWindow = mainWindow
        QWidget.__init__(self)
        self.layout = QGridLayout()
        self.layout.setVerticalSpacing(10)
        self.layout.setHorizontalSpacing(5)
        self.layout.setContentsMargins(0, 10, 60, 10)

        self.previousButton = CustomButton("Previous Run", mainWindow)
        self.nextButton = CustomButton("Next Run", mainWindow)
        self.removePeak = CustomButton("Disable Peak", mainWindow, 1)
        self.updateSigData = CustomButton("Update Sig Vars", mainWindow, 1)
        self.optimizeButton = CustomButton("Fit Sigmoid", mainWindow,1)
        self.calKappaButton = CustomButton("Calculate Kappa", mainWindow, 1)
        self.addSecond = CustomButton("-1 second", mainWindow)
        self.subSecond = CustomButton("+1 second", mainWindow)

        self.removePeak.clicked.connect(self.IgnorePeakClicked)
        self.updateSigData.clicked.connect(self.updateSigVarsClicked)
        self.previousButton.clicked.connect(self.previousButtonClicked)
        self.nextButton.clicked.connect(self.nextButtonClicked)
        self.optimizeButton.clicked.connect(self.optimizeButtonClicked)
        self.calKappaButton.clicked.connect(self.calKappaButtonClicked)
        self.addSecond.clicked.connect(self.addSecondClicked)
        self.subSecond.clicked.connect(self.subSecondClicked)

        self.layout.setColumnStretch(0,0)
        self.layout.addWidget(self.removePeak, 2, 1, 3, 1, alignment = 2)
        self.layout.addWidget(self.updateSigData, 2, 2, 3, 1, alignment = 1)
        self.layout.addWidget(self.previousButton,4,3,4,1)
        self.layout.addWidget(self.nextButton,4,4,4,1)
        self.layout.addWidget(self.subSecond, 0, 3, 3, 1)
        self.layout.addWidget(self.addSecond, 0, 4, 3, 1)
        self.layout.addWidget(self.optimizeButton,2,5,3,1, alignment = 2)
        self.layout.addWidget(self.calKappaButton, 2, 6, 3, 1, alignment = 0)
        self.layout.addWidget(QWidget(),2,7,3,1)

        for i in range(self.layout.columnCount()):
            self.layout.setColumnStretch(i,1)
            self.layout.setColumnMinimumWidth(i,self.width() / 8)
        for i in range(self.layout.rowCount()):
            self.layout.setRowStretch(i, 1)
            self.layout.setRowMinimumHeight(i,self.height()/7)
        self.setLayout(self.layout)

        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Background, settings.controlAreaBackgroundColor)
        self.setPalette(palette)

    def IgnorePeakClicked(self):
        self.mainWindow.controller.disablePeak()

    def updateSigVarsClicked(self):
        updateDialog = InputForm(self.mainWindow)
        if updateDialog.exec_() == QDialog.Accepted:
            (a,b,c) = updateDialog.getData()
            self.mainWindow.controller.reOptimization(a,b,c)

    def addSecondClicked(self):
        self.mainWindow.controller.shiftOneSecond()

    def subSecondClicked(self):
        self.mainWindow.controller.shiftOneSecond(forward=False)

    def nextButtonClicked(self):
        currPeak = self.mainWindow.controller.currPeak
        maxPeak = self.mainWindow.controller.maxPeak
        self.mainWindow.controller.switchToPeak(min(currPeak + 1, maxPeak - 1))

    def previousButtonClicked(self):
        currPeak = self.mainWindow.controller.currPeak
        self.mainWindow.controller.switchToPeak(max(0, currPeak - 1))

    def optimizeButtonClicked(self):
        self.mainWindow.controller.optimizationProcedure()

    def calKappaButtonClicked(self):
        controller = self.mainWindow.controller
        kappaVars = (controller.sigma, controller.temp, controller.dd, controller.iKappa, controller.dd2, controller.iKappa2,controller.solubility)
        sig = kappaVars[0]
        temp = kappaVars[1]
        dd1 = kappaVars[2]
        iKappa1 = kappaVars[3]
        dd2 = kappaVars[4]
        iKappa2 = kappaVars[5]
        solu = kappaVars[6]
        confirmVarDialog = KappaVarConfirmDialog(sig,temp,dd1,iKappa1,dd2,iKappa2,solu,self.mainWindow)
        if confirmVarDialog.exec_() == QDialog.Accepted:
            (sig,temp,dd1,iKappa1,dd2,iKappa2,solu) = confirmVarDialog.getData()
            self.mainWindow.updateKappaVars(sig,temp,dd1,iKappa1,dd2,iKappa2,solu)
            self.mainWindow.calKappa()


class PeakDpDnLogGraphWidget(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.dpView = PeakSquareFigureCanvas()
        self.dNlogView = PeakSquareFigureCanvas()
        self.layout = QHBoxLayout()
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.addWidget(self.dpView)
        self.layout.addWidget(self.dNlogView)
        self.setLayout(self.layout)

    def resize(self, parentWidth, parentHeight):
        self.setFixedWidth(parentWidth)
        self.setFixedHeight(parentHeight * 1 / 2)
        self.dpView.resize(self.width(), self.height())
        self.dNlogView.resize(self.width(), self.height())

class PeakRectFigureCanvas(FigureCanvas):
    def resize(self, parentWidth, parentHeight):
        self.setFixedHeight(parentHeight * 2 / 5 + 3)
        self.setFixedWidth(parentWidth)

    def __init__(self, mainWindow=None):
        fig = Figure(facecolor=settings.graphBackgroundColor)
        super(self.__class__, self).__init__(fig)

    def updateFigure(self, figure):
        self.figure = figure
        self.draw()

        # A hack to make the figure update to the size of the Figure Canvas
        h = self.height()
        self.setFixedHeight(h / 2)
        self.setFixedHeight(h)


class PeakSquareFigureCanvas(FigureCanvas):
    def resize(self, parentWidth, parentHeight):
        self.setFixedHeight(parentHeight)
        self.setFixedWidth(parentWidth / 2)

    def __init__(self, mainWindow=None):
        fig = Figure(facecolor=settings.graphBackgroundColor)
        super(self.__class__, self).__init__(fig)

    def updateFigure(self, figure):
        self.figure = figure
        self.draw()
        h = self.height()
        self.setFixedHeight(h / 2)
        self.setFixedHeight(h)


