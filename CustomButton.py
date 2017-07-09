import sys
import os
from PySide.QtGui import *
from PySide.QtCore import *
import settings
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import time


class CustomButton(QPushButton):
    def __init__(self, text,parent=None,style=0):
        self.parent = parent
        QPushButton.__init__(self, text)
        font = QFont()
        self.setFont(font)
        self.setFlat(True)
        self.setAutoFillBackground(True)
        palette = self.palette()
        if style == 1:
            palette.setColor(QPalette.Button, settings.specialButtonColor)
            palette.setColor(QPalette.ButtonText, settings.simpleButtonTextColor)
        else:
            palette.setColor(QPalette.Button, settings.simpleButtonColor)
            palette.setColor(QPalette.ButtonText, settings.simpleButtonTextColor)
        self.setPalette(palette)

    def resize(self, parent_width, parent_height):
        self.setFixedHeight(parent_height / 3)
        self.setFixedWidth(parent_width / 9)
        font = QFont()
        size = max(1,self.height() * 2 / 11)
        font.setPointSize(size)
        self.setFont(font)
