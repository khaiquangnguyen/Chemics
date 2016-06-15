
import sys
import os
import tempfile

from PySide.QtGui import *
from PySide.QtCore import *
from settings import *

qt_app = QApplication(sys.argv)

width  = 0
height = 0


class MainWindow(QMainWindow):

    def __init__(self,controller):
        QMainWindow.__init__(self)

        #declare a bunch of necessary variables

        #the peak to show graph
        self.peak = 0
        #the progress of the progam
        self.progress = 0
        #the controller of the view
        self.controller = controller
        #the progress dialog
        self.progressDialog = None
        #the gobal width and height for resize of all items inside
        os.chdir(self.controller.tempFolder)
        print os.getcwd()
        global width, height



        # set Size
        self.setWindowTitle('Chemics')
        self.setMinimumHeight(100)
        self.setMinimumWidth(100)
        width = self.width()
        height = self.height()

        #add Menu
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')

        #add select folder button
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
        self.setCentralWidget(ControlPanel())

    def folderSelection(self):
        folder = QFileDialog.getExistingDirectory()
        self.controller.setFolder(folder)


    def showProgress(self):
        if self.progressDialog is None:
            self.progressDialog = QProgressDialog("Tasks in progress...", "Cancel", 0, 10, self)
        self.progressDialog.setWindowModality(Qt.WindowModal)
        self.progressDialog.show()

    def makeProgress(self, message = None):
        self.progress += 1
        self.progressDialog.setValue(self.progress)
        if message is not None:
            self.progressDialog.setLabelText(message)


    def run(self):
        # Show the form
        self.show()
        # Run the qt application
        qt_app.exec_()

    def resizeEvent(self, resizeEvent):
        self.centralWidget().resize()

    def setPeak(self,peak):
        self.peak = peak

    def printMessage(self,message):
        self.centralWidget().subWidget.infoArea.addMessage(message)


class messageArea(QTableWidget):
    def resize(self, height):
        self.setFixedHeight(height * 2 / 5)

    def __init__(self):
        QTableWidget.__init__(self)
        self.setColumnCount(1)
        self.setRowCount(0)
        self.verticalHeader().setVisible(False)
        self.horizontalHeader().setVisible(False)
        self.horizontalHeader().setStretchLastSection(True)
        self.setWordWrap(True)
        self.resizeRowsToContents()

    def addMessage(self, message):
        ##### add code here to process the message before printing
        item = QTableWidgetItem(message)
        self.insertRow(self.rowCount())
        self.setItem(self.rowCount()-1,0,item)

class controlArea(QWidget):
    def resize(self, height):
        self.setFixedHeight(height * 3 / 5)

    def __init__(self):
        QWidget.__init__(self)
        self.layout = QVBoxLayout()
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Background, Qt.red)
        self.setPalette(palette)

class ControlPanel(QWidget):
    ''' An example of PySide/PyQt absolute positioning; the main window
        inherits from QWidget, a convenient widget for an empty window. '''

    def __init__(self):
        # Initialize the object as a QWidget and
        # set its title and minimum width
        QWidget.__init__(self)

        #remove border

        # addSubWidget
        self.subWidget = infoWidget()
        self.graphWidget = graphWidget()
        self.layout = QHBoxLayout()
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.addWidget(self.subWidget)
        self.layout.addWidget(self.graphWidget)
        self.setLayout(self.layout)

    def resize(self):
        self.subWidget.resize(self.width(),self.height())
        self.graphWidget.resize(self.width(),self.height())


class infoWidget(QWidget):
    def __init__(self, parent=None):
        super(self.__class__,self).__init__(parent)
        self.layout = QVBoxLayout()  #Vertical layout
        self.setLayout(self.layout)

        #Set color for layout for easier managing
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Background, Qt.red)
        self.setPalette(palette)
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.controlArea = controlArea()
        self.infoArea = messageArea()
        self.layout.addWidget(self.controlArea)
        self.layout.addWidget(self.infoArea)

    def resize(self,width, height):
        self.setFixedWidth(width / 4 )
        self.setFixedHeight(height)
        self.controlArea.resize(self.height())
        self.infoArea.resize(self.height())


class graphWidget(QWidget):
    def resize(self, width, height):
        self.setFixedWidth(width * 3 / 4)
        self.setFixedHeight(height)
        self.totalView.resize(self.height())

    def __init__(self, parent=None):
        super(self.__class__, self).__init__(parent)
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        # Set color for layout for easier managing
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Background, Qt.blue)
        self.setPalette(palette)

        graph = QPixmap("graph.png")

        self.totalView = totalView()
        self.individualViews = individualViews()
        self.layout.addWidget(self.individualViews)
        self.layout.addWidget(self.totalView)


class totalView(QGraphicsView):
    def resize(self, height):
        self.setFixedHeight(height * 2/5)

    def __init__(self, parent=None):
        super(self.__class__, self).__init__(parent)
        graph = QPixmap(fullGraphFileName)
        self.totalScene = QGraphicsScene()
        self.totalScene.addPixmap(graph)
        self.setScene(self.totalScene)

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


class dpView(QGraphicsView):
    def resize(self, width):
        self.setFixedWidth(width /2 )

    def __init__(self, parent=None):
        super(self.__class__, self).__init__(parent)
        graph = QPixmap(dpViewFileName)
        self.totalScene = QGraphicsScene()
        self.totalScene.addPixmap(graph)
        self.setScene(self.totalScene)
        self.setRenderHint(QPainter.Antialiasing or QPainter.SmoothPixmapTransform or QPainter.HighQualityAntialiasing)
        self.fitInView(self.totalScene.itemsBoundingRect(), Qt.KeepAspectRatio)


class dNlogView(QGraphicsView):
    def resize(self, width):
        self.setFixedWidth(width / 2)

    def __init__(self, parent=None):
        super(self.__class__, self).__init__(parent)
        graph = QPixmap(dnLogFileName)
        self.totalScene = QGraphicsScene()
        self.totalScene.addPixmap(graph)
        self.setScene(self.totalScene)
        self.setRenderHint(QPainter.Antialiasing or QPainter.SmoothPixmapTransform or QPainter.HighQualityAntialiasing)

