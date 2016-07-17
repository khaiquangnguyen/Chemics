
import sys
import os
from PySide.QtGui import *
from PySide.QtCore import *
import settings
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import time

from PeakAlignUI import *
from KappaUI import *

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
        self.mainWindow = mainWindow
        self.infoWidget = InfoWidget(self.mainWindow)
        self.graphWidget = graphWidget(self.mainWindow)
        self.layout = QHBoxLayout()
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.addWidget(self.infoWidget)
        self.layout.addWidget(self.graphWidget)
        self.setLayout(self.layout)

    def resize(self):
        self.infoWidget.resize(self.width(), self.height())
        self.graphWidget.resize(self.width(),self.height())

    def switchToKappa(self):
        self.infoWidget = KappaInfoWidget(self.mainWindow)
        self.graphWidget = KappaGraphWidget(self.mainWindow)
        self.clearLayout(self.layout)
        self.layout.addWidget(self.infoWidget)
        self.layout.addWidget(self.graphWidget)
        self.setLayout(self.layout)
        self.mainWindow.controller.makeKappaGraph()
        self.graphWidget.totalView.updateFigure(self.mainWindow.controller.kappaGraph)

        self.resize()

    def switchToPeak(self):
        self.clearLayout(self.layout)
        self.infoWidget = InfoWidget(self.mainWindow)
        self.graphWidget = graphWidget(self.mainWindow)
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



