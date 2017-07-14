###############################
#
# IMPORT STATEMENTS
#
###############################

import sys
import os
from PySide.QtGui import *
from PySide.QtCore import *
import settings
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import time
from AlignmentUI import *
from KappaUI import *
import webbrowser
import urllib2
from timeit import default_timer as timer

matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4'] = 'PySide'
qt_app = QApplication(sys.argv)

###############################
#
# VARIABLES
#
###############################

width = 0
height = 0
VERSION = 101


###############################
#
# MAIN VIEW CLASS
#
###############################


class View(QMainWindow):
    def __init__(self, controller):
        global VERSION
        global width, height
        QMainWindow.__init__(self)
        self.progress = 0
        self.controller = controller
        self.progress_dialog = None
        self.setWindowTitle('Chemics')
        self.setMinimumHeight(800)
        self.setMinimumWidth(800)
        width = self.width()
        height = self.height()

        # Add navigation menu
        navigation_bar = self.menuBar()
        file_button = navigation_bar.addMenu('&File')
        select_button = QAction('&Select Files', self)
        select_button.setShortcut('Ctrl+O')
        select_button.setStatusTip('Select the files which contains the processed_data files')
        select_button.triggered.connect(self.select_files)
        file_button.addAction(select_button)
        feedback_button = QAction('&Feedback', self)
        feedback_button.triggered.connect(self.submit_feedback)
        navigation_bar.addAction(feedback_button)
        exit_button = QAction('&Exit', self)
        exit_button.setShortcut('Ctrl+Q')
        exit_button.setStatusTip('Exit application')
        exit_button.triggered.connect(self.close)
        file_button.addAction(exit_button)
        self.setCentralWidget(ControlPanel(self))
        self.showMaximized()

    def show_ui(self):
        """
        show the UI of the program.
        """
        self.show()
        qt_app.exec_()

    def show_update_dialog(self, update=1):
        """
        show the update
        """
        updateDialog = QMessageBox()
        if update == 1:
            updateDialog.setText("The program has an update. Please download the update for the program.")
        else:
            updateDialog.setText("The program is up-to-experiment_date.")
        updateDialog.exec_()

    def submit_feedback(self):
        """
        Submit feedback by showing a google form
        """
        webbrowser.open("https://goo.gl/forms/X9OB6AQSJSiKScBs2")

    def check_for_update(self):
        response = urllib2.urlopen('http://khaiquangnguyen.github.io/chemics_update.html')
        html = response.read()
        update = int(html[0])
        self.show_update_dialog(update)

    def select_files(self):
        """
        Select files with processed_data.
        """
        dialog = QFileDialog()
        files = dialog.getOpenFileNames()[0]
        if files:
            self.controller.files = files
            self.controller.run()

    def move_progress_bar_forward(self, message=None, max_value=None, complete=False, value=1):
        """
        Activate the progress bar
        :param max_value: the maximum of the progress dialog. Used to signal reset progress dialog
        :param complete: whether the progress is the completion of the whole procedure
        :param message: the message to show on the progress dialog
        """

        if max_value is not None:
            self.progress = 0
            if self.progress_dialog is not None:
                self.progress_dialog.reset()
                self.progress_dialog.setValue(0)
                self.progress_dialog.setRange(0, max_value)
                if message is not None:
                    self.progress_dialog.setLabelText(message)
                self.progress_dialog.setWindowModality(Qt.WindowModal)
                self.progress_dialog.show()
            else:
                self.progress_dialog = QProgressDialog("Tasks in progress...", "Cancel", 0, max_value, self)
                self.progress_dialog.canceled.connect(self.cancel_progress_bar)
                self.progress_dialog.setWindowModality(Qt.WindowModal)
                self.progress_dialog.show()

        else:
            if self.progress_dialog is not None:
                if complete is True:
                    self.progress_dialog.setValue(self.progress_dialog.maximum())
                    self.progress_dialog.reset()
                else:
                    if message is not None:
                        self.progress_dialog.setLabelText(message)
                    self.progress += 1
                    self.progress_dialog.setValue(self.progress)
                qApp.processEvents()

    def show_error_dialog(self, error_message='Unknown Error!'):
        """
        Show the error message
        :param error_message: The message to show in the error message
        """
        if self.progress_dialog is not None:
            self.progress_dialog.reset()
            self.progress_dialog = None
        warning = QMessageBox()
        warning.setIcon(QMessageBox.Warning)
        warning.setText(error_message)
        warning.exec_()

    def cancel_progress_bar(self):
        """
        Action when the cancelling_progress_bar button of the progress bar is clicked
        """
        self.controller.cancel_progress_bar()

    def resizeEvent(self, resizeEvent):
        self.centralWidget().resize()

    def update_alignment_and_sigmoid_fit_figures(self, scan_time_figure, sigmoid_figure):
        self.centralWidget().graph_widget.alignment_and_sigmoid_fit_view.alignment_view.update_figure(scan_time_figure)
        self.centralWidget().graph_widget.alignment_and_sigmoid_fit_view.sigmoid_fit_view.update_figure(sigmoid_figure)

    def update_temp_and_min_figure(self, new_figure):
        self.centralWidget().graph_widget.temp_and_alignment_view.update_figure(new_figure)

    def calculate_all_kappa_values(self):
        self.centralWidget().switch_to_kappa_widget()
        self.controller.calculate_all_kappa_values()
        self.controller.update_kappa_info_and_graph()

    def update_experiment_information(self):
        self.centralWidget().info_widget.information_table.update_experiment_information()

    def update_scan_information(self):
        self.centralWidget().info_widget.information_table.update_scan_information()
        curr_scan = self.controller.current_scan
        if self.controller.min_pos_SMPS_list[curr_scan] is None or self.controller.min_pos_CCNC_list[curr_scan] is None:
            self.centralWidget().graph_widget.buttons_widget.change_scan_status_button.setText("Enable")
        else:
            if (self.controller.usable_for_sigmoid_fit_list[curr_scan]):
                self.centralWidget().graph_widget.buttons_widget.change_scan_status_button.setText("Disable Scan")
            else:
                self.centralWidget().graph_widget.buttons_widget.change_scan_status_button.setText("Enable Scan")

    def update_scan_information_after_sigmoid_fit(self):
        self.centralWidget().info_widget.information_table.update_scan_information_after_sigmoid_fit()
        curr_scan = self.controller.current_scan
        if self.controller.min_pos_SMPS_list[curr_scan] is None or self.controller.min_pos_CCNC_list[curr_scan] is None:
            self.centralWidget().graph_widget.buttons_widget.change_scan_status_button.setText("Enable Scan")
        else:
            if self.controller.usable_for_kappa_cal_list[curr_scan] and self.controller.usable_for_sigmoid_fit_list[
                                                                                                            curr_scan]:
                self.centralWidget().graph_widget.buttons_widget.change_scan_status_button.setText("Disable Scan")
            else:
                self.centralWidget().graph_widget.buttons_widget.change_scan_status_button.setText("Enable Scan")

    def reset(self):
        self.centralWidget().switch_to_scan_widget()

    def update_kappa_values(self, sigma, temp, dd1, i1, dd2, i2, solu):
        self.controller.sigma = sigma
        self.controller.temp = temp
        self.controller.dd = dd1
        self.controller.iKappa = i1
        self.controller.dd2 = dd2
        self.controller.iKappa2 = i2
        self.controller.solubility = solu

    def update_kappa_info_and_graph(self):
        self.centralWidget().graph_widget.update_figure(self.controller.kappa_graph)
        self.centralWidget().info_widget.update_data()
        self.centralWidget().resize()
        self.centralWidget().info_widget.information_table.toggle_color(self.controller.current_point)
        ss = self.controller.kappa_points_data_list[self.controller.current_point][1]
        dp = self.controller.kappa_points_data_list[self.controller.current_point][0]
        if self.controller.is_show_all_k_points:
            status = self.controller.kappa_points_is_included_list[(dp,ss)]
            if status:
                self.centralWidget().graph_widget.control_widget.toggle_k_point_status_button.setText("Disable")
            else:
                self.centralWidget().graph_widget.control_widget.toggle_k_point_status_button.setText("Enable")

    def get_concentration(self):
        return_value = self.controller.flow_rate
        while True:
            input = QInputDialog.getDouble(self, self.tr("Get Flow Rate"), self.tr("Flow Rate (L/min)"), 0.3,
                                           decimals=5)
            if input[1] is True:
                return_value = float(input[0])
                break
        return return_value


#######################################################
#
# CONTROL PANEL CLASS. CONTAINS ALL THE CONTROL BUTTONS OF THE UI
#
#######################################################

class ControlPanel(QWidget):
    def __init__(self, main_window=None):
        QWidget.__init__(self)
        # addSubWidget
        self.main_window = main_window
        self.layout = QHBoxLayout()
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.info_widget = ScanInformationWidget(self.main_window)
        self.graph_widget = ScanGraphsWidget(self.main_window)
        self.layout.addWidget(self.info_widget)
        self.layout.addWidget(self.graph_widget)
        self.setLayout(self.layout)

    def resize(self):
        self.info_widget.resize(self.width(), self.height())
        self.graph_widget.resize(self.width(), self.height())

    def switch_to_kappa_widget(self):
        self.clear_layout(self.layout)
        self.info_widget = KappaInformationAndDataWidget(self.main_window)
        self.graph_widget = KappaGraphWidget(self.main_window)
        self.layout.addWidget(self.info_widget)
        self.layout.addWidget(self.graph_widget)
        self.setLayout(self.layout)
        self.resize()

    def switch_to_scan_widget(self):
        self.clear_layout(self.layout)
        self.info_widget = ScanInformationWidget(self.main_window)
        self.graph_widget = ScanGraphsWidget(self.main_window)
        self.layout.addWidget(self.info_widget)
        self.layout.addWidget(self.graph_widget)
        self.setLayout(self.layout)
        self.resize()

    def clear_layout(self, layout):
        if layout is not None:
            while layout.count():
                item = layout.takeAt(0)
                widget = item.widget()
                if widget is not None:
                    widget.deleteLater()
                else:
                    self.clear_layout(item.layout())

    def keyReleaseEvent(self, event):
        print event.key()
        if not self.main_window.controller.kappa_ax:
            return
        else:
            self.main_window.controller.on_key_release_kappa_graph(event)



