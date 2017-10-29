<<<<<<< HEAD
from InfoTableStyle import *
from CustomButton import *
import sys
import os
from PySide.QtGui import *
from PySide.QtCore import *
import settings
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import time
from InputDialog import *
from KappaVarConfirmDialog import *

class PeakTextDataWidget(QWidget):
    def __init__(self, mainWindow=None):
        super(self.__class__,self).__init__(mainWindow)
        self.layout = QVBoxLayout()  #Vertical layout
        self.setLayout(self.layout)
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.infoTable = PeakDataTable(mainWindow)
        self.layout.addWidget(self.infoTable)

    def resize(self, parentWidth, parentHeight):
        self.setFixedWidth(parentWidth / 4)
        self.setFixedHeight(parentHeight)
        self.infoTable.resize(self.width(), self.height())

class PeakDataTable(QTableWidget):
    def resize(self, parentWidth, parentHeight):
        self.setFixedHeight(parentHeight)
        self.setFixedWidth(parentWidth)
        self.verticalHeader().setDefaultSectionSize(self.height() / 25)

    def __init__(self, mainWindow = None):
        QTableWidget.__init__(self)
        self.setColumnCount(1)
        self.setRowCount(0)
        self.setShowGrid(False)
        self.verticalHeader().setVisible(False)
        self.horizontalHeader().setVisible(False)
        self.horizontalHeader().setStretchLastSection(True)
        self.verticalHeader().setDefaultSectionSize(self.height() / 25)
        self.setWordWrap(True)
        self.resizeRowsToContents()
        self.setFrameStyle(QFrame.NoFrame)
        self.mainWindow = mainWindow

        #set background color
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, settings.infoAreaBackgroundColor)
        self.setPalette(palette)

    def updateGeneralInfo(self):
        header = TableHeader("Data Information")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        self.addMessage("Date",self.mainWindow.controller.date)
        self.addMessage("Time frame",self.mainWindow.controller.start_time_list[0] + " to " + self.mainWindow.controller.end_time_list[-1])
        self.addMessage("Time per run",self.mainWindow.controller.scan_time)
        self.addMessage("Total run",self.mainWindow.controller.number_of_peak)
        self.addMessage("CPC", self.mainWindow.controller.flowRate)

    def updateBasicPeakInfo(self):
        if self.rowCount() > 6:
            for i in range(self.rowCount() - 6):
                self.removeRow(self.rowCount()-1)
        header = TableHeader("Basic Peak Information")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        currPeak = self.mainWindow.controller.currPeak
        if self.mainWindow.controller.minPosCCNCList[currPeak] and self.mainWindow.controller.minPosCCNCList[currPeak]:
            self.addMessage("Status", "Valid")
        else:
            self.addMessage("Status", "Invalid", color='#EF5350')
        self.addMessage("Current run",currPeak + 1)
        self.addMessage("Saturation",self.mainWindow.controller.superSaturation)

    def updateSigFitPeakInfo(self):
        if self.rowCount() > 6:
            for i in range(self.rowCount() - 6):
                self.removeRow(self.rowCount() - 1)

        # Add basic information data
        header = TableHeader("Single Run Information")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        currPeak = self.mainWindow.controller.currPeak
        if self.mainWindow.controller.minPosCCNCList[currPeak] and self.mainWindow.controller.minPosCCNCList[currPeak]:
            self.addMessage("Status", "Valid")
        else:
            self.addMessage("Status", "Invalid",  color='#EF5350')
        self.addMessage("Current run", currPeak + 1)
        self.addMessage("Saturation", self.mainWindow.controller.superSaturation)

        # Add advance information data
        header = TableHeader("Advance Peak Information")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        if self.mainWindow.controller.usableForKappaCalList[currPeak]:
            self.addMessage("Status", "Valid for Kappa Cal")
        else:
            self.addMessage("Status", "Invalid for Kappa Cal",  color='#EF5350')
        self.addMessage('minDp',self.mainWindow.controller.minDp)
        self.addMessage('minDpAsym', self.mainWindow.controller.minDpAsym)
        self.addMessage('maxDpAsym', self.mainWindow.controller.maxDpAsym)
        self.addMessage("dp50", self.mainWindow.controller.dp50)
        self.addMessage("<dp50 counts",self.mainWindow.controller.dp50LessCount)
        self.addMessage(">dp50 counts",self.mainWindow.controller.dp50MoreCount)
        self.addMessage("dp50(Wet)",self.mainWindow.controller.dp50Wet)
        self.addMessage("dp50+20",self.mainWindow.controller.dp50Plus20)


    def addMessage(self, field,message, color = None):
        ##### add code here to process the message before printing
        message = str(message)
        item = TableItem(field, message, color)
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount()-1,0,item)

class PeakGraphWidget(QWidget):
    def resize(self, parentWidth, parentHeight):
        self.setFixedWidth(parentWidth * 3 / 4)
        self.setFixedHeight(parentHeight)
        self.tempAndMinView.resize(self.width(), self.height())
        self.controlWidget.resize(self.width(), self.height())
        self.dpAndDnlogView.resize(self.width(), self.height())

    def __init__(self, mainWindow=None):
        super(self.__class__, self).__init__(mainWindow)
        self.mainWindow = mainWindow
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.tempAndMinView = PeakRectFigureCanvas()
        self.dpAndDnlogView = PeakDpDnLogGraphWidget()
        self.controlWidget = PeakControlTabWidget(mainWindow)
        self.layout.addWidget(self.dpAndDnlogView)
        self.layout.addWidget(self.controlWidget)
        self.layout.addWidget(self.tempAndMinView)

class PeakControlTabWidget(QWidget):
    def resize(self, parentWidth, parentHeight):
        self.setFixedHeight(parentHeight * 1 / 10)
        self.setFixedWidth(parentWidth)
        self.nextButton.resize(self.width(),self.height())
        self.previousButton.resize(self.width(),self.height())
        self.removePeak.resize(self.width(), self.height())
        self.optimizeButton.resize(self.width(), self.height())
        self.updateSigData.resize(self.width(), self.height())
        self.calKappaButton.resize(self.width(), self.height())
        self.addSecond.resize(self.width(), self.height())
        self.subSecond.resize(self.width(), self.height())

    def __init__(self, mainWindow = None):
        self.mainWindow = mainWindow
        QWidget.__init__(self)
        self.layout = QGridLayout()
        self.layout.setVerticalSpacing(10)
        self.layout.setHorizontalSpacing(5)
        self.layout.setContentsMargins(0, 10, 60, 10)

        self.previousButton = CustomButton("Previous Run", mainWindow)
        self.nextButton = CustomButton("Next Run", mainWindow)
        self.removePeak = CustomButton("Disable Run", mainWindow, 1)
        self.updateSigData = CustomButton("Update Sig Vars", mainWindow, 1)
        self.optimizeButton = CustomButton("Fit Sigmoid", mainWindow,1)
        self.calKappaButton = CustomButton("Calculate Kappa", mainWindow, 1)
        self.addSecond = CustomButton("-1 second", mainWindow)
        self.subSecond = CustomButton("+1 second", mainWindow)

        self.removePeak.clicked.connect(self.IgnorePeakClicked)
        self.updateSigData.clicked.connect(self.updateSigVarsClicked)
        self.previousButton.clicked.connect(self.previousButtonClicked)
        self.nextButton.clicked.connect(self.nextButtonClicked)
        self.optimizeButton.clicked.connect(self.optimizeButtonClicked)
        self.calKappaButton.clicked.connect(self.calKappaButtonClicked)
        self.addSecond.clicked.connect(self.addSecondClicked)
        self.subSecond.clicked.connect(self.subSecondClicked)

        self.layout.setColumnStretch(0,0)
        self.layout.addWidget(self.removePeak, 2, 1, 3, 1, alignment = 2)
        self.layout.addWidget(self.updateSigData, 2, 2, 3, 1, alignment = 1)
        self.layout.addWidget(self.previousButton,4,3,4,1)
        self.layout.addWidget(self.nextButton,4,4,4,1)
        self.layout.addWidget(self.subSecond, 0, 3, 3, 1)
        self.layout.addWidget(self.addSecond, 0, 4, 3, 1)
        self.layout.addWidget(self.optimizeButton,2,5,3,1, alignment = 2)
        self.layout.addWidget(self.calKappaButton, 2, 6, 3, 1, alignment = 0)
        self.layout.addWidget(QWidget(),2,7,3,1)

        for i in range(self.layout.columnCount()):
            self.layout.setColumnStretch(i,1)
            self.layout.setColumnMinimumWidth(i,self.width() / 8)
        for i in range(self.layout.rowCount()):
            self.layout.setRowStretch(i, 1)
            self.layout.setRowMinimumHeight(i,self.height()/7)
        self.setLayout(self.layout)

        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Background, settings.controlAreaBackgroundColor)
        self.setPalette(palette)

    def IgnorePeakClicked(self):
        self.mainWindow.controller.disablePeak()

    def updateSigVarsClicked(self):
        updateDialog = InputForm(self.mainWindow)
        if updateDialog.exec_() == QDialog.Accepted:
            (a,b,c) = updateDialog.getData()
            self.mainWindow.controller.reOptimization(a,b,c)

    def addSecondClicked(self):
        self.mainWindow.controller.shiftOneSecond()

    def subSecondClicked(self):
        self.mainWindow.controller.shiftOneSecond(forward=False)

    def nextButtonClicked(self):
        currPeak = self.mainWindow.controller.currPeak
        number_of_peak = self.mainWindow.controller.number_of_peak
        self.mainWindow.controller.switchToPeak(min(currPeak + 1, number_of_peak - 1))

    def previousButtonClicked(self):
        currPeak = self.mainWindow.controller.currPeak
        self.mainWindow.controller.switchToPeak(max(0, currPeak - 1))

    def optimizeButtonClicked(self):
        self.mainWindow.controller.optimizationProcedure()

    def calKappaButtonClicked(self):
        controller = self.mainWindow.controller
        kappaVars = (controller.sigma, controller.temp, controller.dd, controller.iKappa, controller.dd2, controller.iKappa2,controller.solubility)
        sig = kappaVars[0]
        temp = kappaVars[1]
        dd1 = kappaVars[2]
        iKappa1 = kappaVars[3]
        dd2 = kappaVars[4]
        iKappa2 = kappaVars[5]
        solu = kappaVars[6]
        confirmVarDialog = KappaVarConfirmDialog(sig,temp,dd1,iKappa1,dd2,iKappa2,solu,self.mainWindow)
        if confirmVarDialog.exec_() == QDialog.Accepted:
            (sig,temp,dd1,iKappa1,dd2,iKappa2,solu) = confirmVarDialog.getData()
            self.mainWindow.updateKappaVars(sig,temp,dd1,iKappa1,dd2,iKappa2,solu)
            self.mainWindow.calKappa()


class PeakDpDnLogGraphWidget(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.dpView = PeakSquareFigureCanvas()
        self.dNlogView = PeakSquareFigureCanvas()
        self.layout = QHBoxLayout()
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.addWidget(self.dpView)
        self.layout.addWidget(self.dNlogView)
        self.setLayout(self.layout)

    def resize(self, parentWidth, parentHeight):
        self.setFixedWidth(parentWidth)
        self.setFixedHeight(parentHeight * 1 / 2)
        self.dpView.resize(self.width(), self.height())
        self.dNlogView.resize(self.width(), self.height())

class PeakRectFigureCanvas(FigureCanvas):
    def resize(self, parentWidth, parentHeight):
        self.setFixedHeight(parentHeight * 2 / 5 + 3)
        self.setFixedWidth(parentWidth)

    def __init__(self, mainWindow=None):
        fig = Figure(facecolor=settings.graphBackgroundColor)
        super(self.__class__, self).__init__(fig)

    def updateFigure(self, figure):
        self.figure = figure
        self.draw()

        # A hack to make the figure update to the size of the Figure Canvas
        h = self.height()
        self.setFixedHeight(h / 2)
        self.setFixedHeight(h)


class PeakSquareFigureCanvas(FigureCanvas):
    def resize(self, parentWidth, parentHeight):
        self.setFixedHeight(parentHeight)
        self.setFixedWidth(parentWidth / 2)

    def __init__(self, mainWindow=None):
        fig = Figure(facecolor=settings.graphBackgroundColor)
        super(self.__class__, self).__init__(fig)

    def updateFigure(self, figure):
        self.figure = figure
        self.draw()
        h = self.height()
        self.setFixedHeight(h / 2)
        self.setFixedHeight(h)


=======
from InfoTableStyle import *
from CustomButton import *
import sys
import os
from PySide.QtGui import *
from PySide.QtCore import *
import settings
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import time
from InputDialog import *


class ScanInformationWidget(QWidget):
    def __init__(self, main_window=None):
        super(self.__class__, self).__init__(main_window)
        self.layout = QVBoxLayout()
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.information_table = ScanInformationTable(main_window)
        self.layout.addWidget(self.information_table)
        self.setLayout(self.layout)

    def resize(self, parent_width, parent_height):
        self.setFixedWidth(parent_width / 4)
        self.setFixedHeight(parent_height)
        self.information_table.resize(self.width(), self.height())


class ScanInformationTable(QTableWidget):
    def __init__(self, main_window=None):
        QTableWidget.__init__(self)
        self.setColumnCount(1)
        self.setRowCount(0)
        self.setShowGrid(False)
        self.verticalHeader().setVisible(False)
        self.horizontalHeader().setVisible(False)
        self.horizontalHeader().setStretchLastSection(True)
        self.verticalHeader().setDefaultSectionSize(self.height() / 25)
        self.setWordWrap(True)
        self.resizeRowsToContents()
        self.setFrameStyle(QFrame.NoFrame)
        self.main_window = main_window

        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, settings.infoAreaBackgroundColor)
        self.setPalette(palette)

    def resize(self, parent_width, parent_height):
        self.setFixedHeight(parent_height)
        self.setFixedWidth(parent_width)
        self.verticalHeader().setDefaultSectionSize(self.height() / 25)

    def update_experiment_information(self):
        # add header for table
        header = AllTableHeader("Experiment Information")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        # add table cells
        self.add_message("Experiment Date", self.main_window.controller.experiment_date)
        self.add_message("Experiment Start Time", self.main_window.controller.scan_start_time_list[0])
        self.add_message("Time per scan (s)", self.main_window.controller.scan_duration)
        self.add_message("Total number of scan", self.main_window.controller.number_of_scan)
        self.add_message("Flow rate (L/min)", self.main_window.controller.flow_rate)

    def update_scan_information(self):
        # clear the cells to redraw the table. 6 is the number of cells of exp table
        if self.rowCount() > 6:
            for i in range(self.rowCount() - 6):
                self.removeRow(self.rowCount() - 1)
        # add header
        header = AllTableHeader("Scan Information")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        # add table cells
        current_scan = self.main_window.controller.current_scan
        if self.main_window.controller.is_usable_for_sigmoid_fit_list[current_scan]:
            self.add_message("Usability for Sigmoid Fit", "Positive")
        else:
            self.add_message("Usability for Sigmoid Fit", "Negative", color=settings.NEGATIVE_USABILITY_COLOR)
        self.add_message("Scan Start Time (h/m/s)", self.main_window.controller.scan_start_time_list[current_scan])
        self.add_message("Scan #", current_scan + 1)
        self.add_message("CCNC Data Shift (s)", self.main_window.controller.shift_factor_list[current_scan])
        self.add_message("Super Saturation (%)", self.main_window.controller.super_saturation_rate)

    def update_scan_information_after_sigmoid_fit(self):
        if self.rowCount() > 6:
            for i in range(self.rowCount() - 6):
                self.removeRow(self.rowCount() - 1)
        # add header
        header = AllTableHeader("Scan Information")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        current_scan = self.main_window.controller.current_scan
        # add table cell
        if self.main_window.controller.is_usable_for_sigmoid_fit_list[current_scan]:
            self.add_message("Usability for Sigmoid Fit", "Positive")
        else:
            self.add_message("Usability for Sigmoid Fit", "Negative", color=settings.NEGATIVE_USABILITY_COLOR)
        self.add_message("Scan Start Time (h/m/s)", self.main_window.controller.scan_start_time_list[current_scan])
        self.add_message("Scan #", current_scan + 1)
        self.add_message("CCNC Data Shift (s)", self.main_window.controller.shift_factor_list[current_scan])
        self.add_message("Super Saturation (%)", self.main_window.controller.super_saturation_rate)
        # add header
        header = AllTableHeader("Sigmoid Fit Parameters")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        # add table cells
        if current_scan in self.main_window.controller.unfinished_sigmoid_fit_scans_list:
            self.add_message("Usability for Kappa", "Undecided",color=settings.UNDECIDED_USABILITY_COLOR)
        elif self.main_window.controller.is_usable_for_kappa_cal_list[current_scan] and \
                self.main_window.controller.is_usable_for_sigmoid_fit_list[current_scan]:
            self.add_message("Usability for Kappa", "Positive")
        else:
            self.add_message("Usability for Kappa", "Negative", color=settings.NEGATIVE_USABILITY_COLOR)
        self.add_message('minDp (nm)', self.main_window.controller.min_dp)
        self.add_message('minDpAsym (nm)', self.main_window.controller.min_dp_asym)
        self.add_message('maxDpAsym (nm)', self.main_window.controller.max_dp_asym)
        self.add_message("dp50 (nm)", self.main_window.controller.dp50)
        self.add_message("<dp50 counts", self.main_window.controller.dp50_less_count)
        self.add_message(">dp50 counts", self.main_window.controller.dp50_more_count)
        self.add_message("dp50(Wet) (nm)", self.main_window.controller.dp50_wet)
        self.add_message("dp50+20 (nm)", self.main_window.controller.dp50_plus_20)

    def add_message(self, field, message, color=None):
        if type(message) is not str:
            message = '{0:.2f}'.format(message).rstrip('0').rstrip('.')
        message = str(message)
        item = AlignmentTableItem(field, message, color)
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, item)


class ScanGraphsWidget(QWidget):
    def __init__(self, main_window=None):
        super(self.__class__, self).__init__(main_window)
        self.mainWindow = main_window
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.temp_and_alignment_view = RectFigureCanvas()
        self.alignment_and_sigmoid_fit_view = AlignmentAndSigmoidFitWidget()
        self.buttons_widget = ButtonsWidget(main_window)
        self.layout.addWidget(self.alignment_and_sigmoid_fit_view)
        self.layout.addWidget(self.buttons_widget)
        self.layout.addWidget(self.temp_and_alignment_view)

    def resize(self, parent_width, parent_height):
        self.setFixedWidth(parent_width * 3 / 4)
        self.setFixedHeight(parent_height)
        self.temp_and_alignment_view.resize(self.width(), self.height())
        self.buttons_widget.resize(self.width(), self.height())
        self.alignment_and_sigmoid_fit_view.resize(self.width(), self.height())


class ButtonsWidget(QWidget):
    def __init__(self, main_window=None):
        self.main_window = main_window
        QWidget.__init__(self)
        self.layout = QGridLayout()
        self.layout.setVerticalSpacing(10)
        self.layout.setHorizontalSpacing(10)
        self.layout.setContentsMargins(0, 10, 60, 10)

        self.previous_scan_button = CustomButton("Prev Scan", main_window)
        self.next_scan_button = CustomButton("Next Scan", main_window)
        self.change_scan_status_button = CustomButton("Disable Scan", main_window, 1)
        self.update_sigmoid_fit_parameters = CustomButton("Manual Fit", main_window, 1)
        self.sigmoid_fit_button = CustomButton("Fit Sigmoid", main_window, 1)
        self.calc_kappa_button = CustomButton("Calculate Kappa", main_window, 1)
        self.move_backward_one_second_button = CustomButton("-1 second", main_window)
        self.move_forward_one_second_button = CustomButton("+1 second", main_window)

        self.change_scan_status_button.clicked.connect(self.on_click_change_scan_status)
        self.update_sigmoid_fit_parameters.clicked.connect(self.on_click_update_sigmoid_fit_parameters)
        self.previous_scan_button.clicked.connect(self.on_click_previous_scan)
        self.next_scan_button.clicked.connect(self.on_click_next_scan)
        self.sigmoid_fit_button.clicked.connect(self.on_click_fit_sigmoid_line)
        self.calc_kappa_button.clicked.connect(self.on_click_calculate_kappa)
        self.move_backward_one_second_button.clicked.connect(self.on_click_back_one_second)
        self.move_forward_one_second_button.clicked.connect(self.on_click_forward_one_second)
        self.layout.addWidget(self.change_scan_status_button, 2, 2, 3, 1)
        self.layout.addWidget(self.previous_scan_button, 4, 3, 3, 1)
        self.layout.addWidget(self.next_scan_button, 4, 4, 3, 1)
        self.layout.addWidget(self.move_forward_one_second_button, 0, 3, 3, 1)
        self.layout.addWidget(self.move_backward_one_second_button, 0, 4, 3, 1)
        self.layout.addWidget(self.sigmoid_fit_button, 2, 5, 3, 1)

        for i in range(self.layout.rowCount()):
            self.layout.setRowStretch(i, 0)
            self.layout.setRowMinimumHeight(i, self.height() / 7)
        self.layout.setColumnStretch(0, 1)
        self.layout.setColumnStretch(1, 0)
        self.layout.setColumnStretch(2, 0)
        self.layout.setColumnStretch(3, 0)
        self.layout.setColumnStretch(4, 0)
        self.layout.setColumnStretch(5, 0)
        self.layout.setColumnStretch(6, 0)
        self.layout.setColumnStretch(7, 1)
        self.layout.setColumnMinimumWidth(0, self.width() / 8)
        self.layout.setColumnMinimumWidth(1, self.width() / 8)
        self.layout.setColumnMinimumWidth(2, self.width() / 8)
        self.layout.setColumnMinimumWidth(3, self.width() / 8)
        self.layout.setColumnMinimumWidth(4, self.width() / 8)
        self.layout.setColumnMinimumWidth(5, self.width() / 8)
        self.layout.setColumnMinimumWidth(6, self.width() / 8)
        self.layout.setColumnMinimumWidth(7, self.width() / 8)

        self.setLayout(self.layout)
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Background, settings.controlAreaBackgroundColor)
        self.setPalette(palette)

    def resize(self, parent_width, parent_height):
        self.setFixedHeight(parent_height * 1 / 10)
        self.setFixedWidth(parent_width)
        self.next_scan_button.resize(self.width(), self.height())
        self.previous_scan_button.resize(self.width(), self.height())
        self.change_scan_status_button.resize(self.width(), self.height())
        self.sigmoid_fit_button.resize(self.width(), self.height())
        self.update_sigmoid_fit_parameters.resize(self.width(), self.height())
        self.calc_kappa_button.resize(self.width(), self.height())
        self.move_backward_one_second_button.resize(self.width(), self.height())
        self.move_forward_one_second_button.resize(self.width(), self.height())

    def on_click_change_scan_status(self):
        self.main_window.controller.change_scan_status()

    def on_click_update_sigmoid_fit_parameters(self):
        updateDialog = InputForm(self.main_window)
        if updateDialog.exec_() == QDialog.Accepted:
            (a, b, c) = updateDialog.getData()
            self.main_window.controller.refitting_sigmoid_line(a, b, c)

    def on_click_back_one_second(self):
        self.main_window.controller.shift_data_by_one_second()

    def on_click_forward_one_second(self):
        self.main_window.controller.shift_data_by_one_second(forward=False)

    def on_click_next_scan(self):
        current_scan = self.main_window.controller.current_scan
        number_of_scan = self.main_window.controller.number_of_scan
        self.main_window.controller.switch_to_scan(min(current_scan + 1, number_of_scan - 1))

    def on_click_previous_scan(self):
        current_scan = self.main_window.controller.current_scan
        self.main_window.controller.switch_to_scan(max(0, current_scan - 1))

    def on_click_fit_sigmoid_line(self):
        self.main_window.controller.correct_charges_and_fit_sigmoid_all_scans()
        if self.main_window.controller.finish_scan_alignment_and_auto_sig_fit:
            self.sigmoid_fit_button.hide()
            self.layout.addWidget(self.update_sigmoid_fit_parameters, 2, 1, 3, 1)
            self.layout.addWidget(self.calc_kappa_button, 2, 5, 3, 1)

    def on_click_calculate_kappa(self):
        controller = self.main_window.controller
        if len(controller.unfinished_sigmoid_fit_scans_list) > 0:
            controller.view.show_error_dialog("There are scans with no sigmoid fit. Please go back and finish them manually!")
            return
        kappaVars = (
            controller.sigma, controller.temp, controller.dd, controller.iKappa, controller.dd2, controller.iKappa2,
            controller.solubility)
        sig = kappaVars[0]
        temp = kappaVars[1]
        dd1 = kappaVars[2]
        iKappa1 = kappaVars[3]
        dd2 = kappaVars[4]
        iKappa2 = kappaVars[5]
        solu = kappaVars[6]
        confirmVarDialog = KappaVarConfirmDialog(sig, temp, dd1, iKappa1, dd2, iKappa2, solu, self.main_window)
        if confirmVarDialog.exec_() == QDialog.Accepted:
            (sig, temp, dd1, iKappa1, dd2, iKappa2, solu) = confirmVarDialog.getData()
            self.main_window.update_kappa_values(sig, temp, dd1, iKappa1, dd2, iKappa2, solu)
            self.main_window.calculate_all_kappa_values()


class AlignmentAndSigmoidFitWidget(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.alignment_view = SquareFigureCanvas()
        self.sigmoid_fit_view = SquareFigureCanvas()
        self.layout = QHBoxLayout()
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.addWidget(self.alignment_view)
        self.layout.addWidget(self.sigmoid_fit_view)
        self.setLayout(self.layout)

    def resize(self, parent_width, parent_height):
        self.setFixedWidth(parent_width)
        self.setFixedHeight(parent_height * 1 / 2)
        self.alignment_view.resize(self.width(), self.height())
        self.sigmoid_fit_view.resize(self.width(), self.height())


class RectFigureCanvas(FigureCanvas):
    def __init__(self, main_window=None):
        fig = Figure(facecolor=settings.GRAPH_BACKGROUND_COLOR)
        super(self.__class__, self).__init__(fig)

    def resize(self, parent_width, parent_height):
        self.setFixedHeight(parent_height * 2 / 5 + 3)
        self.setFixedWidth(parent_width)

    def update_figure(self, new_figure):
        if new_figure is None:
            return
        if self.figure != new_figure:
            self.figure = new_figure
            h = self.height()
            self.setFixedHeight(h / 2)
            self.setFixedHeight(h)
            self.draw()
            self.flush_events()

        else:
            self.draw()
            self.flush_events()


class SquareFigureCanvas(FigureCanvas):
    def __init__(self, main_window=None):
        fig = Figure(facecolor=settings.GRAPH_BACKGROUND_COLOR)
        super(self.__class__, self).__init__(fig)

    def resize(self, parent_width, parent_height):
        self.setFixedHeight(parent_height)
        self.setFixedWidth(parent_width / 2)

    def update_figure(self, new_figure):
        if new_figure is None:
            return
        if self.figure != new_figure:
            self.figure = new_figure
            h = self.height()
            self.setFixedHeight(h / 2)
            self.setFixedHeight(h)
            self.draw()
            self.flush_events()

        else:
            self.draw()
            self.flush_events()


class KappaVarConfirmDialog(QDialog):
    def __init__(self, sigma, temp, dd1, i_kappa_1, dd2, i_kappa_2, solu, main_window=None):
        super(self.__class__, self).__init__()
        self.mainWindow = main_window
        self.formLayout = QFormLayout()
        self.sigmaLine = QLineEdit(str(sigma))
        self.tempLine = QLineEdit(str(temp))
        self.dd1Line = QLineEdit(str(dd1))
        self.iKappa1Line = QLineEdit(str(i_kappa_1))
        self.dd2Line = QLineEdit(str(dd2))
        self.iKappa2Line = QLineEdit(str(i_kappa_2))
        self.soluLine = QLineEdit(str(solu))
        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.formLayout.addRow(self.tr("&Sigma"), self.sigmaLine)
        self.formLayout.addRow(self.tr("&Temperature"), self.tempLine)
        self.formLayout.addRow(self.tr("&dry diameter(1)"), self.dd1Line)
        self.formLayout.addRow(self.tr("&iKappa(1)"), self.iKappa1Line)
        self.formLayout.addRow(self.tr("&dry diameter(2)"), self.dd2Line)
        self.formLayout.addRow(self.tr("&iKappa(2)"), self.iKappa2Line)
        self.formLayout.addRow(self.tr("&solubility"), self.soluLine)

        self.buttonBox.accepted.connect(self.acceptInput)
        self.buttonBox.rejected.connect(self.reject)

        self.formLayout.addWidget(self.buttonBox)
        self.setLayout(self.formLayout)

    def acceptInput(self):
        try:
            self.sigma = float(self.sigmaLine.text())
            self.temp = float(self.tempLine.text())
            self.dd1 = float(self.dd1Line.text())
            self.iKappa1 = float(self.iKappa1Line.text())
            self.dd2 = float(self.dd2Line.text())
            self.iKappa2 = float(self.iKappa2Line.text())
            self.solu = float(self.soluLine.text())
            self.accept()
        except:
            self.mainWindow.show_error_dialog("Input processed_data not valid.Please input again!")

    def getData(self):
        return self.sigma, self.temp, self.dd1, self.dd2, self.iKappa1, self.iKappa2, self.solu
>>>>>>> ce5b242c238e196dfe20cf2038e62c260c5c95bd
