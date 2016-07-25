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

class PeakAlignDataWidget(QWidget):
    def __init__(self, mainWindow=None):
        super(self.__class__,self).__init__(mainWindow)
        self.layout = QVBoxLayout()  #Vertical layout
        self.setLayout(self.layout)
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.infoArea = PeakAlignDataTable(mainWindow)
        self.layout.addWidget(self.infoArea)

    def resize(self, parentWidth, parentHeight):
        self.setFixedWidth(parentWidth / 4)
        self.setFixedHeight(parentHeight)
        self.infoArea.resize(self.width(),self.height())

class PeakAlignDataTable(QTableWidget):
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

    def updateData(self):
        header = TableHeader("Data Information")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        self.addMessage("Date",self.mainWindow.controller.date)
        self.addMessage("Time frame",self.mainWindow.controller.startTimeEntries[0] + " to " + self.mainWindow.controller.endTimeEntries[-1])
        self.addMessage("Time per run",self.mainWindow.controller.timeFrame)
        self.addMessage("Total run",self.mainWindow.getMaxPeak())

    def updateBasicPeakInfo(self):
        if self.rowCount() > 5:
            for i in range(4):
                self.removeRow(self.rowCount()-1)
        header = TableHeader("Basic Peak Information")
        currPeak = self.mainWindow.getPeak()
        if self.mainWindow.controller.minPosCCNCList[currPeak] and self.mainWindow.controller.minPosCCNCList[currPeak]:
            self.addMessage("Status", "Valid for curve fit")
        else:
            self.addMessage("Status", "Invalid for curve fit")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        self.addMessage("Current run",self.mainWindow.getPeak() + 1)
        self.addMessage("Saturation",self.mainWindow.controller.superSaturation)

    def updateSigFitPeakInfo(self):
        if self.rowCount() > 8:
            for i in range(10):
                self.removeRow(self.rowCount() - 1)
        header = TableHeader("Advance Peak Information")
        currPeak = self.mainWindow.getPeak()
        if self.mainWindow.controller.usablePeakList[currPeak]:
            self.addMessage("Status", "Valid for Kappa Cal")
        else:
            self.addMessage("Status", "Invalid for Kappa Cal")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        self.addMessage('minDp',self.mainWindow.controller.minDp)
        self.addMessage('minDpAsym', self.mainWindow.controller.minDpAsym)
        self.addMessage('maxDp/maxDpAsym', self.mainWindow.controller.maxDpAsym)
        self.addMessage("dp50", self.mainWindow.controller.dp50)
        self.addMessage("<dp50 counts",self.mainWindow.controller.dp50LessCount)
        self.addMessage(">dp50 counts",self.mainWindow.controller.dp50MoreCount)
        self.addMessage("dp50(Wet)",self.mainWindow.controller.dp50Wet)
        self.addMessage("dp50+20",self.mainWindow.controller.dp50Plus20)


    def addMessage(self, field,message):
        ##### add code here to process the message before printing
        message = str(message)
        item = TableItem(field, message)
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount()-1,0,item)

class graphWidget(QWidget):
    def resize(self, parentWidth, parentHeight):
        self.setFixedWidth(parentWidth * 3 / 4)
        self.setFixedHeight(parentHeight)
        self.totalView.resize(self.width(),self.height())
        self.controlArea.resize(self.width(),self.height())
        self.individualViews.resize(self.width(),self.height())

    def __init__(self, mainWindow=None):
        super(self.__class__, self).__init__(mainWindow)
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.totalView = totalView()
        self.individualViews = individualViews()
        self.controlArea = controlArea(mainWindow)
        self.layout.addWidget(self.individualViews)
        self.layout.addWidget(self.controlArea)
        self.layout.addWidget(self.totalView)

class controlArea(QWidget):
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
        self.removePeak = CustomButton("Ignore Peak", mainWindow, 1)
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
        self.mainWindow.controller.removePeak()

    def updateSigVarsClicked(self):
        updateDialog = InputForm()
        if updateDialog.exec_() == QDialog.Accepted:
            x = updateDialog.getData()
            print x

    def addSecondClicked(self):
        self.mainWindow.addSecond()

    def subSecondClicked(self):
        self.mainWindow.subSecond()

    def nextButtonClicked(self):
        self.mainWindow.updatePeak(min(self.mainWindow.getPeak() + 1, self.mainWindow.getMaxPeak() - 1))

    def previousButtonClicked(self):
        self.mainWindow.updatePeak(max(0, self.mainWindow.getPeak() - 1))

    def optimizeButtonClicked(self):
        self.mainWindow.controller.optimizationProcedure()

    def calKappaButtonClicked(self):
        self.mainWindow.centralWidget().switchToKappa()

class totalView(FigureCanvas):
    def resize(self, parentWidth, parentHeight):
        self.setFixedHeight(parentHeight * 2 / 5 + 5)
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


class individualViews(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.dpView = dpView()
        self.dNlogView = dNlogView()
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


class dpView(FigureCanvas):
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


class dNlogView(FigureCanvas):
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

