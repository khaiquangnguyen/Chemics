import sys
import os
from PySide.QtGui import *
from PySide.QtCore import *
import settings
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import time


class AlignmentTableItem(QWidget):
    def __init__(self,field,message, color = None):
        QWidget.__init__(self)
        self.layout = QVBoxLayout()
        margin = self.height() / 400
        self.layout.setContentsMargins(0, margin, 0, margin)
        self.content = AlignmentInnerTableItem(field, message, color)
        self.layout.addWidget(self.content)
        self.setLayout(self.layout)


class AlignmentInnerTableItem(QWidget):
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


class KappaTableItem(QWidget):
    def __init__(self, message):
        QWidget.__init__(self)
        self.layout = QHBoxLayout()
        self.info_text = InfoText(message)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.addWidget(self.info_text)
        self.setLayout(self.layout)
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, settings.NEGATIVE_USABILITY_BUTTON_COLOR)
        self.setPalette(palette)

    def toggle_color(self):
        self.info_text.toggle_color()


class AllTableHeader(QWidget):
    def __init__(self, header):
        QWidget.__init__(self)
        self.layout = QHBoxLayout()
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
        size = max(10, self.height() * 1/5)
        font.setPointSize(size)
        self.setFont(font)

    def resizeEvent(self,event):
        font = self.font()
        font.setStyleStrategy(QFont.PreferAntialias or QFont.PreferQuality)
        size = max(10, self.height() * 1/5)
        font.setPointSize(size)
        self.setFont(font)

    def toggle_color(self):
        pass


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

    def toggle_color(self):
        pass


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

    def toggle_color(self):
        palette = self.palette()
        palette.setColor(QPalette.Base, settings.TABLE_ROW_HIGHLIGHT_COLOR)
        self.setPalette(palette)

