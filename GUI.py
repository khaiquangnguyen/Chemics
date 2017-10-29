
import sys
import os
from PySide.QtGui import *
from PySide.QtCore import *
import settings
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import time
from AlignmentUI import *
from KappaUI import *
import webbrowser
import urllib2

matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'



qt_app = QApplication(sys.argv)
width = 0
height = 0
VERSION = 101

class MainWindow(QMainWindow):

    def __init__(self,controller):
        global VERSION
        QMainWindow.__init__(self)
        self.progress = 0
        self.controller = controller
        self.progressDialog = None
        self.setWindowTitle('Chemics')
        self.setMinimumHeight(800)
        self.setMinimumWidth(800)
        global width, height
        width = self.width()
        height = self.height()

        # Add Menu
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')

        # Add select files button
        folderSelection = QAction('&Select Files', self)
        folderSelection.setShortcut('Ctrl+O')
        folderSelection.setStatusTip('Select the files which contains the data files')
        folderSelection.triggered.connect(self.folderSelection)
        fileMenu.addAction(folderSelection)

        feedbackSelection = QAction('&Feedback', self)
        feedbackSelection.triggered.connect(self.feedbackSelection)
        menubar.addAction(feedbackSelection)


        updateSelection = QAction('&Update', self)
        updateSelection.triggered.connect(self.updateSelection)
        menubar.addAction(updateSelection)

        #add exit button
        exitAction = QAction( '&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.close)
        fileMenu.addAction(exitAction)

        #set central widget
        self.setCentralWidget(ControlPanel(self))
        self.showMaximized()

        #get update
        # response = urllib2.urlopen('http://khaiquangnguyen.github.io/chemics_update.html')
        # html = response.read()
        # update = int(html[0])
        # if update == VERSION:
        #     self.showUpdateDialog()

    def show_UI(self):
        self.show()
        qt_app.exec_()

    def showUpdateDialog(self, update = 1):
        """ show the update"""
        updateDialog = QMessageBox()
        if update == 1:
            updateDialog.setText("The program has an update. Please download the update for the program.")
        else:
            updateDialog.setText("The program is up-to-date.")
        updateDialog.exec_()

    def feedbackSelection(self):
        """
        Submit feedback by showing a google form
        """
        webbrowser.open("https://goo.gl/forms/Cf6YQtdOXAqGx21U2")

    def updateSelection(self):
        response = urllib2.urlopen('http://khaiquangnguyen.github.io/chemics_update.html')
        html = response.read()
        update = int(html[0])
        self.showUpdateDialog(update)

    def folderSelection(self):
        """
        Select the files which store the file
        """
        dialog = QFileDialog()
        files = dialog.getOpenFileNames()[0]
        if files:
            self.controller.files = files
            self.controller.run()

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


    def resizeEvent(self, resizeEvent):
        self.centralWidget().resize()

    def update_dp_dnlog_figures(self, adjustedFigure, diaFigure):
        self.centralWidget().graphWidget.dpAndDnlogView.dpView.updateFigure(adjustedFigure)
        self.centralWidget().graphWidget.dpAndDnlogView.dNlogView.updateFigure(diaFigure)

    def updateTempOrMinFigure(self, aFigure):
        self.centralWidget().graphWidget.tempAndMinView.updateFigure(aFigure)

    def calKappa(self):
        # if self.controller.completedStep >=2:
        self.controller.calKappa()
        self.centralWidget().switchToKappa()
        self.controller.makeKappaGraph()

    def updateGeneralInfo(self):
        self.centralWidget().infoWidget.infoTable.updateGeneralInfo()

    def updateBasicPeakInfo(self):
        self.centralWidget().infoWidget.infoTable.updateBasicPeakInfo()

    def updateSigFitPeakInfo(self):
        self.centralWidget().infoWidget.infoTable.updateSigFitPeakInfo()

    def reset(self):
        self.centralWidget().switchToPeak()

    def updateKappaVars(self,sigma,temp,dd1,i1,dd2,i2,solu):
        self.controller.sigma = sigma
        self.controller.temp = temp
        self.controller.dd = dd1
        self.controller.iKappa = i1
        self.controller.dd2 = dd2
        self.controller.iKappa2 = i2
        self.controller.solubility = solu

    def updateKappaGraph(self):
        self.centralWidget().graphWidget.graphView.updateFigure(self.controller.kappaGraph)

    def InputFlowRate(self):
        returnValue = self.controller.flowRate
        while True:
            input = QInputDialog.getDouble(self, self.tr("Get Flow Rate"),self.tr("Q(flow rate)"),0.3)
            if input[1] == True:
                returnValue = float(input[0])
                break
        return returnValue

class ControlPanel(QWidget):
    def __init__(self, mainWindow = None):
        QWidget.__init__(self)
        # addSubWidget
        self.mainWindow = mainWindow
        self.layout = QHBoxLayout()
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.infoWidget = PeakTextDataWidget(self.mainWindow)
        self.graphWidget = PeakGraphWidget(self.mainWindow)
        self.layout.addWidget(self.infoWidget)
        self.layout.addWidget(self.graphWidget)
        self.setLayout(self.layout)

    def resize(self):
        self.infoWidget.resize(self.width(), self.height())
        self.graphWidget.resize(self.width(),self.height())

    def switchToKappa(self):
        self.clearLayout(self.layout)
        self.infoWidget = KappaTextDataWidget(self.mainWindow)
        self.graphWidget = KappaGraphWidget(self.mainWindow)
        self.clearLayout(self.layout)
        self.layout.addWidget(self.infoWidget)
        self.layout.addWidget(self.graphWidget)
        self.setLayout(self.layout)
        self.infoWidget.updateData()
        self.resize()

    def switchToPeak(self):
        self.clearLayout(self.layout)
        self.infoWidget = PeakTextDataWidget(self.mainWindow)
        self.graphWidget = PeakGraphWidget(self.mainWindow)
        self.layout.addWidget(self.infoWidget)
        self.layout.addWidget(self.graphWidget)
        self.setLayout(self.layout)
        self.resize()

    def clearLayout(self, layout):
        if layout is not None:
            while layout.count():
                item = layout.takeAt(0)
                widget = item.widget()
                if widget is not None:
                    widget.deleteLater()
                else:
                    self.clearLayout(item.layout())



