
import sys
import os
import tempfile

from PySide.QtGui import *
from PySide.QtCore import *
from settings import *
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'



qt_app = QApplication(sys.argv)
width = 0
height = 0


class MainWindow(QMainWindow):

    def __init__(self,controller):
        QMainWindow.__init__(self)
        # The peakIndex to show graph
        self.peakIndex = 0
        # The progress of the program
        self.progress = 0
        # The controller of the view
        self.controller = controller
        # The progress dialog
        self.progressDialog = None

        # Set the temporary working folder for whatever usage
        os.chdir(self.controller.tempFolder)

        # Set size
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
        # dialog.setViewMode(QFileDialog.Detail)
        folder = dialog.getExistingDirectory(options = QFileDialog.Directory)
        if folder:
            print folder
            self.controller.setFolder(folder)
            self.controller.run()

    def showProgress(self):
        """
        Show the progress dialog
        """
        if self.progressDialog is None:
            self.progressDialog = QProgressDialog("Tasks in progress...", "Cancel", 0, 23, self)
        self.progressDialog.setWindowModality(Qt.WindowModal)
        self.progressDialog.show()

    def makeProgress(self, message = None):
        """
        Progress the progress dialog bar
        :param message: the message to show on the progress dialog
        """
        self.progress += 1
        self.progressDialog.setValue(self.progress)
        if message is not None:
            self.progressDialog.setLabelText(message)

    def showData(self):
        self.centralWidget().update(self.peakIndex)

    def showError(self, errorMessage = 'Unknown Error!'):
        self.progressDialog.reset()
        print self.progressDialog
        warning = QMessageBox()
        warning.setIcon(QMessageBox.Warning)
        warning.setText(errorMessage)
        warning.exec_()

    def run(self):
        # Show the form
        self.show()
        # Run the qt application
        qt_app.exec_()

    def resizeEvent(self, resizeEvent):
        """
        When resize
        """
        self.centralWidget().resize()

    def setPeakIndex(self, peak):
        self.peakIndex = peak

    def updateFigures(self,adjustedFigure,diaFigure):
        print "update Figures"
        self.centralWidget().graphWidget.individualViews.dpView.updateFigure(adjustedFigure)
        self.centralWidget().graphWidget.individualViews.dNlogView.updateFigure(diaFigure)

class ControlPanel(QWidget):
    ''' An example of PySide/PyQt absolute positioning; the main window
        inherits from QWidget, a convenient widget for an empty window. '''

    def __init__(self,parent = None):
        QWidget.__init__(self)
        # addSubWidget
        self.subWidget = InfoWidget(parent)
        self.graphWidget = graphWidget(parent)
        self.layout = QHBoxLayout()
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)

        self.layout.addWidget(self.subWidget)
        self.layout.addWidget(self.graphWidget)
        self.setLayout(self.layout)

        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Background, Qt.white)
        self.setPalette(palette)

    def resize(self):
        self.subWidget.resize(self.width(),self.height())
        self.graphWidget.resize(self.width(),self.height())

    def update(self,peakIndex):
        self.graphWidget.updateImage(peakIndex)


class InfoWidget(QWidget):
    def __init__(self, parent=None):
        super(self.__class__,self).__init__(parent)
        self.layout = QVBoxLayout()  #Vertical layout
        self.setLayout(self.layout)
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.infoArea = messageArea()
        self.layout.addWidget(self.infoArea)

    def resize(self, parentWidth, parentHeight):
        self.setFixedWidth(parentWidth / 4)
        self.setFixedHeight(parentHeight)
        self.infoArea.resize(self.width(),self.height())

class messageArea(QTableWidget):
    def resize(self, parentWidth, parentHeight):
        self.setFixedHeight(parentHeight)
        self.setFixedWidth(parentWidth)

    def __init__(self):
        QTableWidget.__init__(self)
        self.setColumnCount(1)
        self.setRowCount(0)
        self.verticalHeader().setVisible(False)
        self.horizontalHeader().setVisible(False)
        self.horizontalHeader().setStretchLastSection(True)
        self.setWordWrap(True)
        self.resizeRowsToContents()
        self.setFrameStyle(QFrame.NoFrame)

        #set background color
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, "#fbfbfb")
        self.setPalette(palette)
        # self.addMessage("121331231321321")

    def addMessage(self, message):
        ##### add code here to process the message before printing
        item = QTableWidgetItem(message)
        self.insertRow(self.rowCount())
        self.setItem(self.rowCount()-1,0,item)


class graphWidget(QWidget):
    def resize(self, parentWidth, parentHeight):
        self.setFixedWidth(parentWidth * 3 / 4)
        self.setFixedHeight(parentHeight)
        self.totalView.resize(self.width(),self.height())
        self.controlArea.resize(self.width(),self.height())
        self.individualViews.resize(self.width(),self.height())

    def __init__(self, parent=None):
        super(self.__class__, self).__init__(parent)
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.totalView = totalView()
        self.individualViews = individualViews()
        self.controlArea = controlArea(parent)
        self.layout.addWidget(self.individualViews)
        self.layout.addWidget(self.controlArea)
        self.layout.addWidget(self.totalView)

    def updateImage(self,peakIndex):
        self.totalView.updateFigure()

class individualViews(QWidget):
    def __init__(self):
        # Initialize the object as a QWidget and
        # set its title and minimum width
        QWidget.__init__(self)

        # remove border
        # addSubWidget
        self.dpView = dpView()
        self.dNlogView = dNlogView()
        self.layout = QHBoxLayout()
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.addWidget(self.dpView)
        self.layout.addWidget(self.dNlogView)
        self.setLayout(self.layout)
        # fig = Figure(figsize=(600, 600), dpi=72, facecolor=(1, 1, 1), edgecolor=(0, 0, 0))
        # ax = fig.add_subplot(111)
        # ax.plot([0, 1])
        # canvas = FigureCanvas(fig)
        # self.layout.addWidget(canvas)

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
#     def __init__(self, text, parent=None):
#         self.parent = parent
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
        self.showMinScale.resize(self.width(), self.height())

    def __init__(self,parent = None):
        self.parent = parent
        QWidget.__init__(self)
        self.layout = QHBoxLayout()
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)

        self.previousButton = CustomButton("Previous Peak", parent)
        self.nextButton = CustomButton("Next Peakm n", parent)
        self.showMinScale = CustomButton("Min Graph ", parent)
        self.previousButton.clicked.connect(self.previousButtonClicked)
        self.nextButton.clicked.connect(self.nextButtonClicked)
        self.layout.addWidget(self.previousButton)
        self.layout.addWidget(self.nextButton)
        self.layout.addWidget(self.showMinScale)
        self.setLayout(self.layout)

        #background color
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Background, "#37474F")
        self.setPalette(palette)

    def nextButtonClicked(self):
        self.parent.peak +=1
        self.parent.peak = min (self.parent.peak, maxPeak-1)
        self.parent.updateFigures()

    def previousButtonClicked(self):
        self.parent.peak -=1
        self.parent.peak = max(0,self.parent.peak)
        self.parent.updateFigures()

class totalView(FigureCanvas):
    def resize(self, parentWidth, parentHeight):
        self.setFixedHeight(parentHeight * 2 / 5)
        self.setFixedWidth(parentWidth)

    def __init__(self, parent=None):
        super(self.__class__, self).__init__(Figure())


    def updateFigure(self, figure):
        self.figure = figure



class dpView(FigureCanvas):
    def resize(self, parentWidth, parentHeight):
        self.setFixedHeight(parentHeight)
        self.setFixedWidth(parentWidth / 2)

    def __init__(self, parent=None):
        super(self.__class__, self).__init__(Figure())

    def updateFigure(self, figure):
        self.figure = figure
        self.draw()


class dNlogView(FigureCanvas):
    def resize(self, parentWidth, parentHeight):
        self.setFixedHeight(parentHeight)
        self.setFixedWidth(parentWidth / 2)

    def __init__(self, parent=None):
        super(self.__class__, self).__init__(Figure())

    def updateFigure(self, figure):
        self.figure = figure
        self.draw()
