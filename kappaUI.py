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


class KappaTextDataWidget(QWidget):
    def __init__(self, mainWindow=None):
        super(self.__class__,self).__init__(mainWindow)
        self.layout = QVBoxLayout()  #Vertical layout
        self.setLayout(self.layout)
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.rawDataTable = KappaRawDataTable(mainWindow)
        self.graphDataTable = KappaGraphDataTable(mainWindow)
        self.layout.addWidget(self.graphDataTable)
        self.showingGraphData = True
        # self.layout.addWidget(self.rawDataTable)


    def resize(self, parentWidth, parentHeight):
        self.setFixedWidth(parentWidth / 4)
        self.setFixedHeight(parentHeight)
        if self.showingGraphData:
            self.graphDataTable.resize(self.width(), self.height())
        else:
            self.rawDataTable.resize(self.width(), self.height())

    def updateData(self):
        if self.showingGraphData:
            self.graphDataTable.updateData()
        else:
            self.rawDataTable.updateData()

    def changeToGraphDataTable(self):
        self.layout.removeWidget(self.rawDataTable)
        self.layout.addWidget(self.graphDataTable)
        self.showingGraphData = True

    def changeToRawDataTable(self):
        self.layout.removeWidget(self.graphDataTable)
        self.layout.addWidget(self.rawDataTable)
        self.showingGraphData = False

class KappaRawDataTable(QTableWidget):
    def resize(self, parentWidth, parentHeight):
        self.setFixedHeight(parentHeight)
        self.setFixedWidth(parentWidth)
        self.verticalHeader().setDefaultSectionSize(self.height() / 25)

    def __init__(self, mainWindow = None):
        QTableWidget.__init__(self)
        self.setColumnCount(5)
        self.setRowCount(0)
        self.setShowGrid(False)
        self.verticalHeader().setVisible(False)
        self.horizontalHeader().setVisible(False)
        self.horizontalHeader().setStretchLastSection(True)
        self.verticalHeader().setDefaultSectionSize(self.height() / 25)
        self.setWordWrap(True)
        self.setFrameStyle(QFrame.NoFrame)
        self.mainWindow = mainWindow
        #set background color
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, settings.infoAreaBackgroundColor)
        self.setPalette(palette)

    def updateData(self):
        self.insertRow(self.rowCount())
        ss = TableHeader('ss')
        self.setCellWidget(self.rowCount() - 1, 0, ss)
        dp = TableHeader('dp')
        self.setCellWidget(self.rowCount() - 1, 1, dp)
        app = TableHeader('K/app')
        self.setCellWidget(self.rowCount() - 1, 2, app)
        ana = TableHeader('K/ana')
        self.setCellWidget(self.rowCount() - 1, 3, ana)
        devi = TableHeader('% devi')
        self.setCellWidget(self.rowCount() - 1, 4, devi)

        kappaDict = self.mainWindow.getKappaDict()
        for aKey in kappaDict.keys():
            if kappaDict[aKey]:
                for aSS in kappaDict[aKey]:
                    ss = aKey
                    dp = aSS[0]
                    app = aSS[1]
                    ana = aSS[2]
                    devi = aSS[3]
                    self.addMessage(ss,dp,app,ana,devi)

    def addMessage(self, ss,dp,app,ana,devi):
        ss = ('% .2f' % ss)
        dp = ('% .4f' % dp)
        app =('% .4f' % app)
        ana = ('% .4f' % ana)
        devi = ('% .4f' % devi)
        self.insertRow(self.rowCount())
        ss = SingleTableItem(ss)
        self.setCellWidget(self.rowCount() - 1, 0, ss)
        dp = SingleTableItem(dp)
        self.setCellWidget(self.rowCount() - 1, 1, dp)
        app = SingleTableItem(app)
        self.setCellWidget(self.rowCount() - 1, 2, app)
        ana = SingleTableItem(ana)
        self.setCellWidget(self.rowCount() - 1, 3, ana)
        devi = SingleTableItem(devi)
        self.setCellWidget(self.rowCount() - 1, 4, devi)

class KappaGraphDataTable(QTableWidget):
    def resize(self, parentWidth, parentHeight):
        self.setFixedHeight(parentHeight)
        self.setFixedWidth(parentWidth)
        self.verticalHeader().setDefaultSectionSize(self.height() / 25)

    def __init__(self, mainWindow=None):
        QTableWidget.__init__(self)
        self.setRowCount(0)
        self.setShowGrid(False)
        self.verticalHeader().setVisible(False)
        self.horizontalHeader().setVisible(False)
        self.horizontalHeader().setStretchLastSection(True)
        self.verticalHeader().setDefaultSectionSize(self.height() / 25)
        self.setWordWrap(True)
        self.setFrameStyle(QFrame.NoFrame)

        self.mainWindow = mainWindow

        # set background color
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, settings.infoAreaBackgroundColor)
        self.setPalette(palette)

    def updateData(self):
        dataDict = self.mainWindow.getAlphaPineneDict()
        self.setColumnCount(len(dataDict.keys()) + 1)
        self.headerList = []
        ssHeader = SingleTableHeaderItem("SS(%)")
        meanDPHeader = SingleTableHeaderItem("MeanDP")
        stdDPheader = SingleTableHeaderItem("StdDP")
        meanAppHeader = SingleTableHeaderItem("Mean K, app")
        stdAppHeader = SingleTableHeaderItem("Std K, app")
        meanAnaHeader = SingleTableHeaderItem("Mean K, ana")
        stdAnaHeader = SingleTableHeaderItem("Std K, ana")
        meanDeviHeader = SingleTableHeaderItem("Mean %Deviation")
        deviMeanHeader = SingleTableHeaderItem("%Deviation of Mean")
        self.headerList.append(ssHeader)
        self.headerList.append(meanDPHeader)
        self.headerList.append(stdDPheader)
        self.headerList.append(meanAppHeader)
        self.headerList.append(stdAppHeader)
        self.headerList.append(meanAnaHeader)
        self.headerList.append(stdAnaHeader)
        self.headerList.append(meanDeviHeader)
        self.headerList.append(deviMeanHeader)

        # Insert Header
        for i in range(0,len(self.headerList)):
            self.insertRow(self.rowCount())
            self.setCellWidget(self.rowCount() - 1, 0, self.headerList[i])

        # Insert data
        count = 1
        for aKey in dataDict.keys():
            aList = [aKey]
            aList.extend(dataDict[aKey])
            self.addMessage(aList,count)
            count += 1


    def addMessage(self, dataList, columnPos):
        for i in range(len(dataList)):
            aCell =  ('% .2f' % dataList[i])
            aCell = SingleTableItem(aCell)
            self.setCellWidget(i, columnPos, aCell)




class KappaGraphDataWidget(QWidget):
    def resize(self, parentWidth, parentHeight):
        self.setFixedWidth(parentWidth * 3 / 4)
        self.setFixedHeight(parentHeight)
        self.totalView.resize(self.width(),self.height())
        self.controlArea.resize(self.width(),self.height())

    def __init__(self, mainWindow=None):
        super(self.__class__, self).__init__(mainWindow)
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.totalView = KappaFigureCanvas()
        self.controlArea = KappaControlTabWidget(mainWindow)
        self.layout.addWidget(self.totalView)
        self.layout.addWidget(self.controlArea)


class KappaControlTabWidget(QWidget):
    def resize(self, parentWidth, parentHeight):
        self.setFixedHeight(parentHeight * 1 / 10)
        self.setFixedWidth(parentWidth)
        self.nextButton.resize(self.width(),self.height())
        self.previousButton.resize(self.width(),self.height())
        self.showMinGraphButton.resize(self.width(), self.height())
        self.optimizeButton.resize(self.width(), self.height())
        self.showTotalGraphButton.resize(self.width(), self.height())
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
        self.showMinGraphButton = CustomButton("Minimum Graph ", mainWindow,1)
        self.showTotalGraphButton = CustomButton("Complete Graph", mainWindow,1)
        self.optimizeButton = CustomButton("Fit Sigmoid", mainWindow,1)
        self.calKappaButton = CustomButton("Calculate Kappa", mainWindow, 1)
        self.addSecond = CustomButton("+1 second", mainWindow)
        self.subSecond = CustomButton("-1 second", mainWindow)
        self.previousButton.clicked.connect(self.previousButtonClicked)
        self.nextButton.clicked.connect(self.nextButtonClicked)
        self.optimizeButton.clicked.connect(self.optimizeButtonClicked)
        self.calKappaButton.clicked.connect(self.calKappaButtonClicked)

        self.layout.setColumnStretch(0,0)

        self.layout.addWidget(self.showMinGraphButton, 2, 1, 3, 1,alignment = 2)
        self.layout.addWidget(self.showTotalGraphButton, 2, 2, 3, 1, alignment = 1)

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

    def nextButtonClicked(self):
        self.mainWindow.updatePeak(min(self.mainWindow.getPeak() + 1, self.mainWindow.getMaxPeak() - 2))

    def previousButtonClicked(self):
        self.mainWindow.updatePeak(max(0, self.mainWindow.getPeak() - 1))

    def optimizeButtonClicked(self):
        self.mainWindow.controller.optimizationProcedure()

    def calKappaButtonClicked(self):
        self.mainWindow.centralWidget().switchToKappa()


class KappaFigureCanvas(FigureCanvas):
    def resize(self, parentWidth, parentHeight):
        self.setFixedHeight(parentHeight * 9 / 10 + 5)
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

