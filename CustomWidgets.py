from PySide.QtGui import *
from PySide.QtCore import *


class TitleHLine(QWidget):
    def __init__(self, title):
        super(self.__class__, self).__init__()
        h_layout = QHBoxLayout()
        h_layout.addWidget(QHLine())
        h_layout.addWidget(QLabel(title))
        h_layout.addWidget(QHLine())
        self.setLayout(h_layout)


class QVLine(QFrame):
    def __init__(self):
        super(QVLine, self).__init__()
        self.setFrameShape(QFrame.VLine)
        self.setFrameShadow(QFrame.Sunken)


class QHLine(QFrame):
    def __init__(self):
        super(QHLine, self).__init__()
        self.setFrameShape(QFrame.HLine)
        self.setFrameShadow(QFrame.Raised)


class LabelDoublelSpinbox(QWidget):
    def __init__(self, title):
        super(self.__class__, self).__init__()
        h_layout = QHBoxLayout()
        label = QLabel(title)
        self.content_box = QDoubleSpinBox()
        self.content_box.setMinimumWidth(100)
        self.content_box.setAlignment(Qt.AlignCenter)
        h_layout.addWidget(label)
        h_layout.addWidget(self.content_box)
        h_layout.setAlignment(Qt.AlignRight)
        h_layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(h_layout)

    def setMaximum(self, max_value):
        self.content_box.setMaximum(max_value)

    def setValue(self, n):
        self.content_box.setValue(n)


class LabelSpinbox(QWidget):
    def __init__(self, title):
        super(self.__class__, self).__init__()
        h_layout = QHBoxLayout()
        label = QLabel(title)
        self.content_box = QSpinBox()
        self.content_box.setMinimumWidth(100)
        self.content_box.setAlignment(Qt.AlignCenter)
        h_layout.addWidget(label)
        h_layout.addWidget(self.content_box)
        h_layout.setAlignment(Qt.AlignRight)
        h_layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(h_layout)

    def setMaximum(self, max_value):
        self.content_box.setMaximum(max_value)

    def setValue(self, n):
        self.content_box.setValue(n)


class ArrowSpinBox(QWidget):
    def __init__(self, forward):
        super(self.__class__, self).__init__()
        h_layout = QHBoxLayout()
        self.content_box = QSpinBox()
        self.content_box.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.content_box.setMinimumWidth(100)
        self.content_box.setAlignment(Qt.AlignCenter)
        self.prev_button = QToolButton()
        self.prev_button.setArrowType(Qt.LeftArrow)
        self.next_button = QToolButton()
        self.next_button.setArrowType(Qt.RightArrow)
        if forward:
            self.prev_button.clicked.connect(self.content_box.stepDown)
            self.next_button.clicked.connect(self.content_box.stepUp)
        else:
            self.prev_button.clicked.connect(self.content_box.stepUp)
            self.next_button.clicked.connect(self.content_box.stepDown)
        h_layout.addWidget(self.prev_button)
        h_layout.addWidget(self.content_box)
        h_layout.addWidget(self.next_button)
        h_layout.setAlignment(Qt.AlignRight)
        h_layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(h_layout)

    def setCallback(self, fn):
        self.content_box.valueChanged.connect(fn)

    def setValue(self, value):
        self.content_box.setValue(value)

    def setRange(self, min, max):
        self.content_box.setRange(min, max)


class LabelSpinBox(QWidget):
    def __init__(self, left_text, right_text, forward):
        super(self.__class__, self).__init__()
        h_layout = QHBoxLayout()
        self.content_box = QSpinBox()
        self.content_box.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.content_box.setAlignment(Qt.AlignCenter)
        self.content_box.setMinimumWidth(60)
        self.left_button = QPushButton(left_text)
        self.left_button.setAutoDefault(False)
        self.right_button = QPushButton(right_text)
        self.right_button.setAutoDefault(False)
        if forward:
            self.left_button.clicked.connect(self.content_box.stepDown)
            self.right_button.clicked.connect(self.content_box.stepUp)
        else:
            self.left_button.clicked.connect(self.content_box.stepUp)
            self.right_button.clicked.connect(self.content_box.stepDown)
        h_layout.addWidget(self.left_button)
        h_layout.addWidget(self.content_box)
        h_layout.addWidget(self.right_button)
        h_layout.setAlignment(Qt.AlignRight)
        h_layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(h_layout)

    def setCallback(self, fn):
        self.content_box.valueChanged.connect(fn)

    def setValue(self, value):
        self.content_box.setValue(value)

    def setRange(self, min, max):
        self.content_box.setRange(min, max)


class LongLabel(QLabel):
    def paintEvent(self, event):
        painter = QPainter(self)
        metrics = QFontMetrics(self.font())
        elided = metrics.elidedText(self.text(), Qt.ElideRight, self.width())
        painter.drawText(self.rect(), self.alignment(), elided)


class KappaTableWidget(QTableWidget):
    def __init__(self, kappa_docker):
        super(self.__class__, self).__init__()
        self.kappa_docker = kappa_docker
        self.setColumnCount(4)
        self.setHorizontalHeaderLabels(["Status", "ss", "dp50", "app-K"])
        self.setSelectionBehavior(QTableWidget.SelectRows)
        self.itemClicked.connect(self.rowClick)
        self.itemChanged.connect(self.toogleRow)
        self.itemDoubleClicked.connect(self.rowClick)
        self.setColumnWidth(1, 60)
        self.setColumnWidth(2, 80)
        self.setColumnWidth(3, 80)

    def addRow(self, ss, dp_50, app_k, status):
        self.insertRow(self.rowCount())
        current_row = self.rowCount() - 1
        status_box = QTableWidgetItem()
        if status == True:
            status_box.setCheckState(Qt.Checked)
            status_box.setText("Included")
        else:
            status_box.setCheckState(Qt.Unchecked)
            status_box.setText("Excluded")
        status_box.setTextAlignment(Qt.AlignCenter)
        ss = QTableWidgetItem(str(ss))
        ss.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
        ss.setTextAlignment(Qt.AlignCenter)
        dp_50 = QTableWidgetItem(str(dp_50))
        dp_50.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
        dp_50.setTextAlignment(Qt.AlignCenter)
        app_k = QTableWidgetItem(str(app_k))
        app_k.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
        app_k.setTextAlignment(Qt.AlignCenter)
        self.setItem(current_row, 0, status_box)
        self.setItem(current_row, 1, ss)
        self.setItem(current_row, 2, dp_50)
        self.setItem(current_row, 3, app_k)
        self.sortByColumn(1, Qt.AscendingOrder)
        self.sortByColumn(2, Qt.AscendingOrder)

    def rowClick(self, item):
        # if item is not the checkbox. This is because sometimes, pyqt
        # will recognize everything as checkable
        if item.column() != 0:
            return
        row = item.row()
        if self.item(row, 1) is None or self.item(row, 2) is None:
            return
        ss = self.item(row, 1).text()
        dp = self.item(row, 2).text()
        # self.kappa_docker.select_k_point(ss, dp)

    def toogleRow(self, item):
        row = item.row()
        # only check for the first column
        if item.column() != 0:
            return
        if self.item(row, 1) is None or self.item(row, 2) is None:
            return
        ss = self.item(row, 1).text()
        dp = self.item(row, 2).text()
        if item.checkState() == Qt.Unchecked:
            self.kappa_docker.toggle_k_points(ss, dp, False)
        else:
            self.kappa_docker.toggle_k_points(ss, dp, True)

    def setStatus(self, ss, dp, status):
        # find the right row first
        for i in range(self.rowCount()):
            if float(self.item(i, 1).text()) == ss and float(self.item(i, 2).text()) == dp:
                # this kappa point is valid
                status_box = self.item(i, 0)
                if status == True:
                    status_box.setCheckState(Qt.Checked)
                    status_box.setText("Included")
                else:
                    status_box.setCheckState(Qt.Unchecked)
                    status_box.setText("Excluded")
