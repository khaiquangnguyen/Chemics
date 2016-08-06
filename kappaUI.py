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
        self.mainWindow = mainWindow
        self.layout = QVBoxLayout()  #Vertical layout
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.currDataTable = KappaConstDataTable(mainWindow)
        self.layout.addWidget(self.currDataTable)
        self.setLayout(self.layout)

    def resize(self, parentWidth, parentHeight):
        self.setFixedWidth(parentWidth / 4)
        self.setFixedHeight(parentHeight)
        self.currDataTable.resize(self.width(), self.height())

    def updateData(self):
        self.currDataTable.updateData()
        self.currDataTable.resize(self.width(), self.height())

    def changeToGraphDataTable(self):
        self.clearLayout(self.layout)
        self.currDataTable = KappaGraphDataTable(self.mainWindow)
        self.layout.addWidget(self.currDataTable)
        self.updateData()

    def changeToRawDataTable(self):
        self.clearLayout(self.layout)
        self.currDataTable = KappaRawDataTable(self.mainWindow)
        self.layout.addWidget(self.currDataTable)
        self.updateData()

    def changeToConstDataTable(self):
        self.clearLayout(self.layout)
        self.currDataTable = KappaConstDataTable(self.mainWindow)
        self.layout.addWidget(self.currDataTable)
        self.updateData()

    def clearLayout(self, layout):
        if layout is not None:
            while layout.count():
                item = layout.takeAt(0)
                widget = item.widget()
                if widget is not None:
                    widget.deleteLater()
                else:
                    self.clearLayout(item.layout())

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
        for i in range(self.rowCount()):
            self.removeRow(self.rowCount()-1)

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
        for i in range(self.rowCount()):
            self.removeRow(self.rowCount() - 1)

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

class KappaConstDataTable(QTableWidget):
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
        for i in range(self.rowCount()):
            self.removeRow(self.rowCount() - 1)

        header = TableHeader("Variables")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        kappaVars = self.mainWindow.getKappaVars()
        s = kappaVars[0]
        t = kappaVars[1]
        d1 = kappaVars[2]
        i1 = kappaVars[3]
        d2 = kappaVars[4]
        i2 = kappaVars[5]
        solu = kappaVars[6]
        self.addMessage("Sigma",s)
        self.addMessage("Temperature",t)
        self.addMessage("dry diamater(1)",d1)
        self.addMessage("iKppa(1)",i1)
        self.addMessage("dry diamater(2)", d2)
        self.addMessage("iKppa(2)", i2)
        self.addMessage("Solubility",solu)

    def addMessage(self, field,message):
        ##### add code here to process the message before printing
        message = str(message)
        item = TableItem(field, message)
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount()-1,0,item)

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
        self.graphDataButton.resize(self.width(),self.height())
        self.rawDataButton.resize(self.width(), self.height())
        self.variableButton.resize(self.width(), self.height())
        self.fullGraphButton.resize(self.width(), self.height())
        self.focusedGraphButton.resize(self.width(), self.height())

    def __init__(self, mainWindow = None):
        self.mainWindow = mainWindow
        QWidget.__init__(self)
        self.layout = QHBoxLayout()
        self.layout.setContentsMargins(0, 10, 60, 10)

        self.rawDataButton = CustomButton("Raw Data", mainWindow)
        self.graphDataButton = CustomButton("Graph Data", mainWindow)
        self.variableButton = CustomButton("Variables", mainWindow)
        self.fullGraphButton = CustomButton("Full Graph", mainWindow)
        self.focusedGraphButton = CustomButton("Focused Graph", mainWindow)
        self.rawDataButton.clicked.connect(self.rawDataButtonClicked)
        self.graphDataButton.clicked.connect(self.graphButtonClicked)
        self.variableButton.clicked.connect(self.variableButtonClicked)
        self.fullGraphButton.clicked.connect(self.fullGraphClicked)
        self.focusedGraphButton.clicked.connect(self.focusedGraphClicked)

        self.layout.addWidget(self.rawDataButton)
        self.layout.addWidget(self.graphDataButton)
        self.layout.addWidget(self.variableButton)
        self.layout.addWidget(self.focusedGraphButton)
        self.layout.addWidget(self.fullGraphButton)
        self.setLayout(self.layout)

        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Background, settings.controlAreaBackgroundColor)
        self.setPalette(palette)

    def graphButtonClicked(self):
        self.mainWindow.centralWidget().infoWidget.changeToGraphDataTable()


    def rawDataButtonClicked(self):
        self.mainWindow.centralWidget().infoWidget.changeToRawDataTable()

    def variableButtonClicked(self):
        self.mainWindow.centralWidget().infoWidget.changeToConstDataTable()

    def fullGraphClicked(self):
        self.mainWindow.controller.makeKappaGraph(fullGraph=True)

    def focusedGraphClicked(self):
        self.mainWindow.controller.makeKappaGraph(fullGraph=False)

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

