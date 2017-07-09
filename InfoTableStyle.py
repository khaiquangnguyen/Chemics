import sys
import os
from PySide.QtGui import *
from PySide.QtCore import *
import settings
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import time


class TableItem(QWidget):
    def __init__(self,field,message, color = None):
        QWidget.__init__(self)
        self.layout = QVBoxLayout()
        margin = self.height() / 400
        self.layout.setContentsMargins(0, margin, 0, margin)
        self.content = InnerTableItem(field, message, color)
        self.layout.addWidget(self.content)
        self.setLayout(self.layout)


class InnerTableItem(QWidget):
    def __init__(self,field,message, color):
        QWidget.__init__(self)
        self.layout = QGridLayout()
        self.layout.setHorizontalSpacing(5)
        self.fieldText = FieldText(field)
        self.infoText = InfoText(message, color)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.addWidget(self.fieldText, 0, 0, 1, 6)
        self.layout.addWidget(self.infoText, 0, 7, 1, 4)
        self.setLayout(self.layout)

        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, settings.infoAreaItemBackgroundColor)
        self.setPalette(palette)


class SingleTableItem(QWidget):
    def __init__(self, message):
        QWidget.__init__(self)
        self.layout = QHBoxLayout()
        self.infoText = InfoText(message)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.addWidget(self.infoText)
        self.setLayout(self.layout)
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, settings.infoAreaItemBackgroundColor)
        self.setPalette(palette)

class SingleTableHeaderItem(QWidget):
    def __init__(self, message):
        QWidget.__init__(self)
        self.layout = QHBoxLayout()
        self.fieldText = FieldText(message)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.addWidget(self.fieldText)
        self.setLayout(self.layout)
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, settings.infoAreaItemBackgroundColor)
        self.setPalette(palette)

class TableHeader(QWidget):
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
        size = max(10, self.height() * 2/9)
        font.setPointSize(size)
        self.setFont(font)

    def resizeEvent(self,event):
        font = self.font()
        font.setStyleStrategy(QFont.PreferAntialias or QFont.PreferQuality)
        size = max(10, self.height() * 2/9)
        font.setPointSize(size)
        self.setFont(font)


class FieldText(QLabel):
    def __init__(self,message):
        QLabel.__init__(self,message)
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, settings.infoAreaFieldColor)
        palette.setColor(QPalette.Text, settings.infoAreaFontColor)
        self.setPalette(palette)
        self.setContentsMargins(5, 0, 0, 0)
        font = QFont()
        size = max(4, self.height() * 2 / 10)
        font.setPointSize(size)
        self.setFont(font)

    def resizeEvent(self, event):
        font = self.font()
        size = max(10, self.height() * 2 / 10)
        font.setPointSize(size)
        self.setFont(font)

class InfoText(QLabel):
    def __init__(self, message, color = None):
        QLabel.__init__(self, message)
        self.setAutoFillBackground(True)
        palette = QPalette()
        if not color:
            palette.setColor(QPalette.Base, settings.infoAreaItemColor)
        else:
            palette.setColor(QPalette.Base, color)
        palette.setColor(QPalette.Text, settings.infoAreaFontColor)
        self.setPalette(palette)
        self.setContentsMargins(10, 0, 0, 0)
        font = QFont()
        size = max(4, self.height() * 2 / 10)
        font.setPointSize(size)
        self.setFont(font)

    def resizeEvent(self, event):
        font = self.font()
        size = max(10, self.height() * 2 / 10)
        font.setPointSize(size)
        self.setFont(font)
