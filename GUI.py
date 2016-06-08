
import sys

from PySide.QtGui import *
from PySide.QtCore import *

qt_app = QApplication(sys.argv)

width  = 0
height = 0


class MainWindow(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)

        # set Size
        self.setWindowTitle('Chemics')
        # screen = QDesktopWidget().screenGeometry()
        # width = screen.width()
        # height = screen.height()
        global width, height
        self.setMinimumHeight(100)
        self.setMinimumWidth(100)
        width = self.width()
        height = self.height()
        #add Menu
        exitAction = QAction( '&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.close)
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(exitAction)
        self.setCentralWidget(ControlPanel())


    def run(self):
        # Show the form
        self.show()
        # Run the qt application
        qt_app.exec_()

    def resizeEvent(self, resizeEvent):
        self.centralWidget().resize()

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

    def resize(self,width, height):
        self.setFixedWidth(width / 4 )
        self.setFixedHeight(height)


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
        self.layout.addWidget(self.totalView)
        self.layout.addWidget(self.individualViews)

class totalView(QGraphicsView):
    def resize(self, height):
        self.setFixedHeight(height * 2/5)

    def __init__(self, parent=None):
        super(self.__class__, self).__init__(parent)
        graph = QPixmap("graph.png")
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
        graph = QPixmap("graph.png")
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
        graph = QPixmap("graph.png")
        self.totalScene = QGraphicsScene()
        self.totalScene.addPixmap(graph)
        self.setScene(self.totalScene)
        self.setRenderHint(QPainter.Antialiasing or QPainter.SmoothPixmapTransform or QPainter.HighQualityAntialiasing)



# Create an instance of the application window and run it
app = MainWindow()
app.run()