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


class ScanInformationWidget(QWidget):
    def __init__(self, main_window=None):
        super(self.__class__,self).__init__(main_window)
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
    def __init__(self, main_window = None):
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
        header = TableHeader("Experiment Information")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        self.add_message("Experiment Date", self.main_window.controller.experiment_date)
        self.add_message("Experiment Start Time", self.main_window.controller.scan_start_time_list[0])
        self.add_message("Time per scan (s)", self.main_window.controller.scan_duration)
        self.add_message("Total number of scan", self.main_window.controller.number_of_scan)
        self.add_message("Concentration Rate", self.main_window.controller.concentration)

    def update_scan_information(self):
        if self.rowCount() > 6:
            for i in range(self.rowCount() - 6):
                self.removeRow(self.rowCount()-1)
        header = TableHeader("Scan Information")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        current_scan = self.main_window.controller.current_scan
        if self.main_window.controller.min_pos_CCNC_list[current_scan] and self.main_window.controller.min_pos_CCNC_list[current_scan]:
            self.add_message("Usability for Sigmoid Fit", "Positive")
        else:
            self.add_message("Usability for Sigmoid Fit", "Negative", color='#EF5350')
        self.add_message("Scan #", current_scan + 1)
        self.add_message("Saturation Rate", self.main_window.controller.super_saturation_rate)

    def update_scan_information_after_sigmoid_fit(self):
        if self.rowCount() > 6:
            for i in range(self.rowCount() - 6):
                self.removeRow(self.rowCount() - 1)

        # Add basic information processed_data
        header = TableHeader("Scan Information")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        current_scan = self.main_window.controller.current_scan
        if self.main_window.controller.min_pos_CCNC_list[current_scan] and self.main_window.controller.min_pos_CCNC_list[current_scan]:
            self.add_message("Usability for Sigmoid Fit", "Positive")
        else:
            self.add_message("Usability for Sigmoid Fit", "Negative", color='#EF5350')
        self.add_message("Scan #", current_scan + 1)
        self.add_message("Saturation Rate", self.main_window.controller.super_saturation_rate)

        header = TableHeader("Sigmoid Fit Parameters")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        if self.main_window.controller.usable_for_kappa_cal_list[current_scan]:
            self.add_message("Usability for Calculating Kappa", "Positive")
        else:
            self.add_message("Usability for Calculating Kappa", "Negative", color='#EF5350')
        self.add_message('minDp', self.main_window.controller.min_dp)
        self.add_message('minDpAsym', self.main_window.controller.min_dp_asym)
        self.add_message('maxDpAsym', self.main_window.controller.max_dp_asym)
        self.add_message("dp50", self.main_window.controller.dp50)
        self.add_message("<dp50 counts", self.main_window.controller.dp50_less_count)
        self.add_message(">dp50 counts", self.main_window.controller.dp50_more_count)
        self.add_message("dp50(Wet)", self.main_window.controller.dp50_wet)
        self.add_message("dp50+20", self.main_window.controller.dp50_plus_20)

    def add_message(self, field, message, color=None):
        if type(message) is not str:
            message = '{0:.2f}'.format(message)
        message = str(message)
        item = TableItem(field, message, color)
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount()-1,0,item)


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
        self.layout.setHorizontalSpacing(5)
        self.layout.setContentsMargins(0, 10, 60, 10)

        self.previous_scan_button = CustomButton("Previous Run", main_window)
        self.next_scan_button = CustomButton("Next Run", main_window)
        self.invalid_scan_button = CustomButton("Disable Scan", main_window, 1)
        self.update_sigmoid_fit_parameters = CustomButton("Update Sig Pars", main_window, 1)
        self.sigmoid_fit_button = CustomButton("Fit Sigmoid", main_window, 1)
        self.calc_kappa_button = CustomButton("Calculate Kappa", main_window, 1)
        self.move_backward_one_second_button = CustomButton("-1 second", main_window)
        self.move_forward_one_second_button = CustomButton("+1 second", main_window)

        self.invalid_scan_button.clicked.connect(self.on_click_disable_scan)
        self.update_sigmoid_fit_parameters.clicked.connect(self.on_click_update_sigmoid_fit_parameters)
        self.previous_scan_button.clicked.connect(self.on_click_previous_scan)
        self.next_scan_button.clicked.connect(self.on_click_next_scan)
        self.sigmoid_fit_button.clicked.connect(self.on_click_fit_sigmoid_line)
        self.calc_kappa_button.clicked.connect(self.on_click_calculate_kappa)
        self.move_backward_one_second_button.clicked.connect(self.on_click_back_one_second)
        self.move_forward_one_second_button.clicked.connect(self.on_click_forward_one_second)

        self.layout.setColumnStretch(0,0)
        self.layout.addWidget(self.invalid_scan_button, 2, 1, 3, 1, alignment = 2)
        self.layout.addWidget(self.update_sigmoid_fit_parameters, 2, 2, 3, 1, alignment = 1)
        self.layout.addWidget(self.previous_scan_button, 4, 3, 4, 1)
        self.layout.addWidget(self.next_scan_button, 4, 4, 4, 1)
        self.layout.addWidget(self.move_forward_one_second_button, 0, 3, 3, 1)
        self.layout.addWidget(self.move_backward_one_second_button, 0, 4, 3, 1)
        self.layout.addWidget(self.sigmoid_fit_button, 2, 5, 3, 1, alignment = 2)
        self.layout.addWidget(self.calc_kappa_button, 2, 6, 3, 1, alignment = 0)
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

    def resize(self, parent_width, parent_height):
        self.setFixedHeight(parent_height * 1 / 10)
        self.setFixedWidth(parent_width)
        self.next_scan_button.resize(self.width(), self.height())
        self.previous_scan_button.resize(self.width(), self.height())
        self.invalid_scan_button.resize(self.width(), self.height())
        self.sigmoid_fit_button.resize(self.width(), self.height())
        self.update_sigmoid_fit_parameters.resize(self.width(), self.height())
        self.calc_kappa_button.resize(self.width(), self.height())
        self.move_backward_one_second_button.resize(self.width(), self.height())
        self.move_forward_one_second_button.resize(self.width(), self.height())

    def on_click_disable_scan(self):
        self.main_window.controller.invalidate_scan()

    def on_click_update_sigmoid_fit_parameters(self):
        updateDialog = InputForm(self.main_window)
        if updateDialog.exec_() == QDialog.Accepted:
            (a,b,c) = updateDialog.getData()
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

    def on_click_calculate_kappa(self):
        controller = self.main_window.controller
        kappaVars = (controller.sigma, controller.temp, controller.dd, controller.iKappa, controller.dd2, controller.iKappa2,controller.solubility)
        sig = kappaVars[0]
        temp = kappaVars[1]
        dd1 = kappaVars[2]
        iKappa1 = kappaVars[3]
        dd2 = kappaVars[4]
        iKappa2 = kappaVars[5]
        solu = kappaVars[6]
        confirmVarDialog = KappaVarConfirmDialog(sig, temp, dd1, iKappa1, dd2, iKappa2, solu, self.main_window)
        if confirmVarDialog.exec_() == QDialog.Accepted:
            (sig,temp,dd1,iKappa1,dd2,iKappa2,solu) = confirmVarDialog.getData()
            self.main_window.update_kappa_values(sig, temp, dd1, iKappa1, dd2, iKappa2, solu)
            self.main_window.calculate_kappa_value()


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
        fig = Figure(facecolor=settings.graphBackgroundColor)
        super(self.__class__, self).__init__(fig)

    def resize(self, parent_width, parent_height):
        self.setFixedHeight(parent_height * 2 / 5 + 3)
        self.setFixedWidth(parent_width)

    def update_figure(self, new_figure):
        self.figure = new_figure
        self.draw()
        h = self.height()
        self.setFixedHeight(h / 2)
        self.setFixedHeight(h)


class SquareFigureCanvas(FigureCanvas):
    def __init__(self, main_window=None):
        fig = Figure(facecolor=settings.graphBackgroundColor)
        super(self.__class__, self).__init__(fig)

    def resize(self, parent_width, parent_height):
        self.setFixedHeight(parent_height)
        self.setFixedWidth(parent_width / 2)

    def update_figure(self, new_figure):
        self.figure = new_figure
        self.draw()
        h = self.height()
        self.setFixedHeight(h / 2)
        self.setFixedHeight(h)


