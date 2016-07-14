
import sys
import os
from PySide.QtGui import *
from PySide.QtCore import *
from settings import *
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import time

matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'



qt_app = QApplication(sys.argv)
width = 0
height = 0


class MainWindow(QMainWindow):

    def __init__(self,controller):
        QMainWindow.__init__(self)
        self.progress = 0
        self.controller = controller
        self.progressDialog = None
        os.chdir(self.controller.tempFolder)
        self.setWindowTitle('Chemics')
        self.setMinimumHeight(800)
        self.setMinimumWidth(800)
        global width, height
        width = self.width()
        height = self.height()

        # Add Menu
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')

        # Add select folder button
        folderSelection = QAction('&Select Folder', self)
        folderSelection.setShortcut('Ctrl+O')
        folderSelection.setStatusTip('Select the folder which contains the data files')
        folderSelection.triggered.connect(self.folderSelection)
        fileMenu.addAction(folderSelection)

        #add exit button
        exitAction = QAction( '&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.close)
        fileMenu.addAction(exitAction)

        #set central widget
        self.setCentralWidget(ControlPanel(self))
        self.showMaximized()

    def folderSelection(self):
        """
        Select the folder which store the file
        """
        dialog = QFileDialog()
        dialog.setFileMode(QFileDialog.Directory)
        folder = dialog.getExistingDirectory(options = QFileDialog.Directory)
        if folder:
            self.controller.setFolder(folder)

    def makeProgress(self, message = None, maxValue = None, complete = False, value = 1):
        """
        Progress the progress dialog bar
        :param maxValue: the maximum of the progress dialog. Used to signal reset progress dialog
        :param complete: whether the progress is the completion of the whole procedure
        :param message: the message to show on the progress dialog
        """

        if maxValue is not None:
            self.progress = 0
            if self.progressDialog is not None:
                self.progressDialog.reset()
                self.progressDialog.setValue(0)
                self.progressDialog.setRange(0,maxValue)
                if message is not None:
                    self.progressDialog.setLabelText(message)
                self.progressDialog.setWindowModality(Qt.WindowModal)
                self.progressDialog.show()
            else:
                self.progressDialog = QProgressDialog("Tasks in progress...", "Cancel", 0, maxValue, self)
                self.progressDialog.canceled.connect(self.cancelProgress)
                self.progressDialog.setWindowModality(Qt.WindowModal)
                self.progressDialog.show()

        else:
            if self.progressDialog is not None:
                if complete == True:
                    self.progressDialog.setValue(self.progressDialog.maximum())
                    self.progressDialog.reset()
                else:
                    if message is not None:
                        self.progressDialog.setLabelText(message)
                    self.progress += 1
                    self.progressDialog.setValue(self.progress)
                qApp.processEvents()

    def showError(self, errorMessage = 'Unknown Error!'):
        """
        Show the error message
        :param errorMessage: The message to show in the error message
        """
        if self.progressDialog is not None:
            self.progressDialog.reset()
            self.progressDialog = None
        warning = QMessageBox()
        warning.setIcon(QMessageBox.Warning)
        warning.setText(errorMessage)
        warning.exec_()

    def cancelProgress(self):
        """
        Action when the cancel button of the progress bar is clicked
        """
        self.controller.cancelProgress()

    def run(self):
        self.show()
        qt_app.exec_()

    def resizeEvent(self, resizeEvent):
        self.centralWidget().resize()

    def updateFigures(self,adjustedFigure,diaFigure):
        """
        Update the two figures in the adjust-figure area and diameter-figure area
        :param adjustedFigure: the figure to update into the area
        :param diaFigure: the figure to update into the area
        """
        self.centralWidget().graphWidget.individualViews.dpView.updateFigure(adjustedFigure)
        self.centralWidget().graphWidget.individualViews.dNlogView.updateFigure(diaFigure)

    def updateGeneralInfo(self):
        """
        Update the information area
        """
        self.centralWidget().infoWidget.infoArea.updateGeneralInfo()

    def updatePeakInfo(self):
        self.centralWidget().infoWidget.infoArea.updatePeakInfo()

    def updatePeak(self,peak):
        self.controller.updatePeak(peak)

    def getPeak(self):
        return self.controller.currPeak

    def getMaxPeak(self):
        return self.controller.maxPeak

    def updateTotalViewFigure(self, aFigure):
        """
        Update the figure in the total view. Either the graph of the whole experiment or of the minimum
        :param aFigure: The figure to update
        """
        self.centralWidget().graphWidget.totalView.updateFigure(aFigure)

class ControlPanel(QWidget):
    def __init__(self, mainWindow = None):
        QWidget.__init__(self)
        # addSubWidget
        self.infoWidget = InfoWidget(mainWindow)
        self.graphWidget = graphWidget(mainWindow)

        self.layout = QHBoxLayout()
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.addWidget(self.infoWidget)
        self.layout.addWidget(self.graphWidget)
        self.setLayout(self.layout)

        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Background, Qt.white)
        self.setPalette(palette)

    def resize(self):
        self.infoWidget.resize(self.width(), self.height())
        self.graphWidget.resize(self.width(),self.height())

class InfoWidget(QWidget):
    def __init__(self, mainWindow=None):
        super(self.__class__,self).__init__(mainWindow)
        self.layout = QVBoxLayout()  #Vertical layout
        self.setLayout(self.layout)
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.infoArea = messageArea(mainWindow)
        self.layout.addWidget(self.infoArea)

    def resize(self, parentWidth, parentHeight):
        self.setFixedWidth(parentWidth / 4)
        self.setFixedHeight(parentHeight)
        self.infoArea.resize(self.width(),self.height())

class messageArea(QTableWidget):
    def resize(self, parentWidth, parentHeight):
        self.setFixedHeight(parentHeight)
        self.setFixedWidth(parentWidth)
        self.verticalHeader().setDefaultSectionSize(self.height() / 15)

    def __init__(self, mainWindow = None):
        QTableWidget.__init__(self)
        self.setColumnCount(1)
        self.setRowCount(0)
        self.setShowGrid(False)
        self.verticalHeader().setVisible(False)
        self.horizontalHeader().setVisible(False)
        self.horizontalHeader().setStretchLastSection(True)
        self.verticalHeader().setDefaultSectionSize(self.height() / 15)
        self.setWordWrap(True)
        self.resizeRowsToContents()
        self.setFrameStyle(QFrame.NoFrame)
        self.mainWindow = mainWindow


        #set background color
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, "#fbfbfb")
        self.setPalette(palette)

    def updateGeneralInfo(self):
        header = sectionHeader("Data Information")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        self.addMessage(self.mainWindow.controller.folder)
        self.addMessage(self.mainWindow.controller.date)
        self.addMessage(self.mainWindow.controller.startTimeEntries[0] + "--> " + self.mainWindow.controller.endTimeEntries[-1])
        self.addMessage(self.mainWindow.controller.timeFrame)
        self.addMessage(self.mainWindow.controller.maxPeak - 1)
        header = sectionHeader("Peak Information")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)

    def updatePeakInfo(self):
        if self.rowCount() > 7:
            for i in range(7):
                self.removeRow(self.rowCount()-1)
        print self.rowCount()
        self.addMessage(self.mainWindow.getPeak())
        self.addMessage(self.mainWindow.controller.superSaturation)
        self.addMessage(self.mainWindow.controller.dp50)
        self.addMessage(self.mainWindow.controller.dp50LessCount)
        self.addMessage(self.mainWindow.controller.dp50MoreCount)
        self.addMessage(self.mainWindow.controller.dp50Wet)
        self.addMessage(self.mainWindow.controller.dp50Plus20)


    def addMessage(self, message):
        ##### add code here to process the message before printing
        message = str(message)
        item = outerTableItem(message)
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount()-1,0,item)


class outerTableItem(QWidget):
    def __init__(self,message):
        QWidget.__init__(self)
        self.layout = QVBoxLayout()
        margin = self.height() / 100
        self.layout.setContentsMargins(0, margin, 0, margin)
        self.insideItem = tableItem(message)
        self.layout.addWidget(self.insideItem)
        self.setLayout(self.layout)

class tableItem(QWidget):
    def __init__(self,message):
        QWidget.__init__(self)
        self.layout = QGridLayout()
        self.layout.setHorizontalSpacing(5)
        self.fieldText = fieldTextCustom(message)
        self.infoText = infoTextCustom(message)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.addWidget(self.fieldText,0,0,1,1)
        self.layout.addWidget(self.infoText,0,1,1,2)
        self.setLayout(self.layout)

        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, "#000000")
        self.setPalette(palette)

class sectionHeader(QWidget):
    def __init__(self, header):
        QWidget.__init__(self)
        self.layout = QVBoxLayout()
        margin = self.height() / 100
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.header = QLabel(header)
        self.layout.addWidget(self.header)
        self.setLayout(self.layout)


class fieldTextCustom(QLabel):
    def __init__(self,message):
        QLabel.__init__(self,message)
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, "#ffffff")
        self.setPalette(palette)

class infoTextCustom(QLabel):
    def __init__(self, message):
        QLabel.__init__(self, message)
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, "#ffffff")
        self.setPalette(palette)


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


class CustomButton(QPushButton):
    def __init__(self, text,parent=None):
        self.parent = parent
        QPushButton.__init__(self, text)
        font = QFont()
        self.setFont(font)
        self.setFlat(True)
        self.setAutoFillBackground(True)

    def resize(self,parentWidth, parentHeight):
        self.setFixedHeight(parentHeight * 3 / 8 )
        self.setFixedWidth(parentWidth / 12)
        font = QFont()
        size = max(5,self.height() / 3.5)
        font.setPointSize(size)
        self.setFont(font)


# class CustomLabel(QLabel):
#     def __init__(self, text, mainWindow=None):
#         self.mainWindow = mainWindow
#         QLabel.__init__(self, text)
#         font = QFont()
#         self.setFont(font)
#
#     def resize(self, parentWidth, parentHeight):
#         self.setFixedHeight(parentHeight * 3 / 8)
#         self.setFixedWidth(parentWidth / 12)
#         font = QFont()
#         size = max(5, self.height() * 3/ 4)
#         font.setPointSize(size)
#         self.setFont(font)
#

class controlArea(QWidget):
    def resize(self, parentWidth, parentHeight):
        self.setFixedHeight(parentHeight * 1 / 10)
        self.setFixedWidth(parentWidth)
        self.nextButton.resize(self.width(),self.height())
        self.previousButton.resize(self.width(),self.height())
        self.showMinGraphButton.resize(self.width(), self.height())

    def __init__(self, mainWindow = None):
        self.mainWindow = mainWindow
        QWidget.__init__(self)
        self.layout = QHBoxLayout()
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)

        self.previousButton = CustomButton("Previous Peak", mainWindow)
        self.nextButton = CustomButton("Next Peak", mainWindow)
        self.showMinGraphButton = CustomButton("Min Graph ", mainWindow)
        self.showTotalGraphButton = CustomButton("Full Graph", mainWindow)
        self.optimizeButton = CustomButton("Optimize", mainWindow)

        self.previousButton.clicked.connect(self.previousButtonClicked)
        self.nextButton.clicked.connect(self.nextButtonClicked)
        self.optimizeButton.clicked.connect(self.optimize)
        self.layout.addWidget(self.previousButton)
        self.layout.addWidget(self.nextButton)
        self.layout.addWidget(self.showMinGraphButton)
        self.layout.addWidget(self.optimizeButton)
        self.setLayout(self.layout)

        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Background, "#37474F")
        self.setPalette(palette)

    def nextButtonClicked(self):
        self.mainWindow.updatePeak(min(self.mainWindow.getPeak() + 1, self.mainWindow.getMaxPeak() - 2))

    def previousButtonClicked(self):
        self.mainWindow.updatePeak(max(0, self.mainWindow.getPeak() - 1))

    def optimize(self):
        self.mainWindow.controller.optimizationProcedure()

class totalView(FigureCanvas):
    def resize(self, parentWidth, parentHeight):
        self.setFixedHeight(parentHeight * 2 / 5)
        self.setFixedWidth(parentWidth)

    def __init__(self, mainWindow=None):
        super(self.__class__, self).__init__(Figure())

    def updateFigure(self, figure):
        self.figure = figure
        self.draw()

        # A hack to make the figure update to the size of the Figure Canvas
        h = self.height()
        self.setFixedHeight(h / 2)
        self.setFixedHeight(h)

class dpView(FigureCanvas):
    def resize(self, parentWidth, parentHeight):
        self.setFixedHeight(parentHeight)
        self.setFixedWidth(parentWidth / 2)

    def __init__(self, mainWindow=None):
        super(self.__class__, self).__init__(Figure())

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
        super(self.__class__, self).__init__(Figure())

    def updateFigure(self, figure):
        self.figure = figure
        self.draw()
        h = self.height()
        self.setFixedHeight(h / 2)
        self.setFixedHeight(h)

