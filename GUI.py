import PySide
import sys

from PySide.QtGui import *
from PySide.QtGui import *
from PySide.QtCore import *

qt_app = QApplication(sys.argv)


class LayoutExample(QWidget):
    ''' An example of PySide/PyQt absolute positioning; the main window
        inherits from QWidget, a convenient widget for an empty window. '''

    def __init__(self):
        # Initialize the object as a QWidget and
        # set its title and minimum width
        QWidget.__init__(self)
        self.setWindowTitle('Chemics')
        screen = QDesktopWidget().screenGeometry()
        self.setMinimumSize(screen.width()/2,screen.height()/2)

        # addSubWidget
        self.subWidget = infoWidget()
        self.graphWidget = graphWidget()
        self.layout = QGridLayout()
        self.layout.addWidget(self.subWidget,0,0,0,1)
        self.layout.addWidget(self.graphWidget,0,1,0,3)
        self.setLayout(self.layout)


    def run(self):
        # Show the form
        self.show()
        # Run the qt application
        qt_app.exec_()

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


class graphWidget(QWidget):
    def __init__(self, parent=None):
        super(self.__class__, self).__init__(parent)
        self.layout = QVBoxLayout()  # Vertical layout
        self.setLayout(self.layout)

        # Set color for layout for easier managing
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Background, Qt.blue)
        self.setPalette(palette)

# Create an instance of the application window and run it
app = LayoutExample()
app.run()