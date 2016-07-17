import sys
import os
from PySide.QtGui import *
from PySide.QtCore import *
import settings
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import time


class outerTableItem(QWidget):
    def __init__(self,field,message):
        QWidget.__init__(self)
        self.layout = QVBoxLayout()
        margin = self.height() / 400
        self.layout.setContentsMargins(0, margin, 0, margin)
        self.insideItem = tableItem(field,message)
        self.layout.addWidget(self.insideItem)
        self.setLayout(self.layout)

class tableItem(QWidget):
    def __init__(self,field,message):
        QWidget.__init__(self)
        self.layout = QGridLayout()
        self.layout.setHorizontalSpacing(5)
        self.fieldText = fieldTextCustom(field)
        self.infoText = infoTextCustom(message)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.addWidget(self.fieldText,0,0,1,1)
        self.layout.addWidget(self.infoText,0,1,1,3)
        self.setLayout(self.layout)

        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, settings.infoAreaItemBackgroundColor)
        self.setPalette(palette)

class sectionHeader(QWidget):
    def __init__(self, header):
        QWidget.__init__(self)
        self.layout = QVBoxLayout()
        self.layout.setContentsMargins(5, 0, 0, 0)
        self.header = QLabel(header)
        palette = self.header.palette()
        palette.setColor(QPalette.Text, settings.infoAreaHeaderFontColor)
        self.header.setPalette(palette)

        self.layout.addWidget(self.header)
        self.setLayout(self.layout)

        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, settings.infoAreaHeaderColor)
        self.setPalette(palette)

        font = QFont()
        size = max(10, self.height() * 3 / 9)
        font.setPointSize(size)
        self.setFont(font)

    def resizeEvent(self,event):
        font = self.font()
        font.setStyleStrategy(QFont.PreferAntialias or QFont.PreferQuality)
        size = max(10, self.height() * 3 / 9)
        font.setPointSize(size)
        self.setFont(font)

class fieldTextCustom(QLabel):
    def __init__(self,message):
        QLabel.__init__(self,message)
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, settings.infoAreaFieldColor)
        palette.setColor(QPalette.Text, settings.infoAreaFontColor)
        self.setPalette(palette)
        self.setContentsMargins(5,0,0,0)
        font = QFont()
        size = max(4, self.height() * 3 / 10)
        font.setPointSize(size)
        self.setFont(font)

    def resizeEvent(self, event):
        font = self.font()
        size = max(10, self.height() * 3 / 10)
        font.setPointSize(size)
        self.setFont(font)

class infoTextCustom(QLabel):
    def __init__(self, message):
        QLabel.__init__(self, message)
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, settings.infoAreaItemColor)
        palette.setColor(QPalette.Text, settings.infoAreaFontColor)
        self.setPalette(palette)
        self.setContentsMargins(10, 0, 0, 0)
        font = QFont()
        size = max(4, self.height() * 3 / 10)
        font.setPointSize(size)
        self.setFont(font)

    def resizeEvent(self, event):
        font = self.font()
        size = max(10, self.height() * 3 / 10)
        font.setPointSize(size)
        self.setFont(font)
