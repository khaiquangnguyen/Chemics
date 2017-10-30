<<<<<<< HEAD
=======
<<<<<<< HEAD
=======
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5
###############################
#
# IMPORT STATEMENTS
#
###############################
<<<<<<< HEAD
=======
>>>>>>> ce5b242c238e196dfe20cf2038e62c260c5c95bd
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5

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
<<<<<<< HEAD
=======
<<<<<<< HEAD

matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'


=======
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5
from timeit import default_timer as timer

matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4'] = 'PySide'
# qt_app = QApplication(sys.argv)

###############################
#
# VARIABLES
#
###############################
<<<<<<< HEAD
=======
>>>>>>> ce5b242c238e196dfe20cf2038e62c260c5c95bd
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5

qt_app = QApplication(sys.argv)
width = 0
height = 0
VERSION = 1


<<<<<<< HEAD
class View(QMainWindow):
    def __init__(self, controller):
=======
class MainWindow(QMainWindow):

<<<<<<< HEAD
    def __init__(self,controller):
=======

class View(QMainWindow):
    def __init__(self, controller):
>>>>>>> ce5b242c238e196dfe20cf2038e62c260c5c95bd
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5
        global VERSION
        QMainWindow.__init__(self)
        self.progress = 0
        self.controller = controller
        self.progress_dialog = None
        self.setWindowTitle('Chemics')
        self.setMinimumHeight(800)
        self.setMinimumWidth(800)
        global width, height
        width = self.width()
        height = self.height()

<<<<<<< HEAD
=======
<<<<<<< HEAD
        # Add Menu
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')

        # Add select files button
        folderSelection = QAction('&Select Files', self)
        folderSelection.setShortcut('Ctrl+O')
        folderSelection.setStatusTip('Select the files which contains the data files')
        folderSelection.triggered.connect(self.folderSelection)
        fileMenu.addAction(folderSelection)

        feedbackSelection = QAction('&Feedback', self)
        feedbackSelection.triggered.connect(self.feedbackSelection)
        menubar.addAction(feedbackSelection)


        updateSelection = QAction('&Update', self)
        updateSelection.triggered.connect(self.updateSelection)
        menubar.addAction(updateSelection)

        #add exit button
        exitAction = QAction( '&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.close)
        fileMenu.addAction(exitAction)

        #set central widget
=======
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5
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
<<<<<<< HEAD
=======
>>>>>>> ce5b242c238e196dfe20cf2038e62c260c5c95bd
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5
        self.setCentralWidget(ControlPanel(self))
        self.showMaximized()

        #get update
        # response = urllib2.urlopen('http://khaiquangnguyen.github.io/chemics_update.html')
        # html = response.read()
        # update = int(html[0])
        # if update == VERSION:
        #     self.showUpdateDialog()

    def show_ui(self):
        self.show()
        # self.check_for_update()
        qt_app.exec_()

<<<<<<< HEAD
=======
<<<<<<< HEAD
    def showUpdateDialog(self, update = 1):
        """ show the update"""
        updateDialog = QMessageBox()
        if update == 1:
            updateDialog.setText("The program has an update. Please download the update for the program.")
        else:
            updateDialog.setText("The program is up-to-date.")
        updateDialog.exec_()
=======
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5
    def show_update_dialog(self, update=1):
        """
        show the update
        """
        if update > VERSION:
            updateDialog = QMessageBox()
            updateDialog.setText("The program has an update. Please download the update for the program.")
            updateDialog.exec_()
<<<<<<< HEAD
=======
>>>>>>> ce5b242c238e196dfe20cf2038e62c260c5c95bd
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5

    def submit_feedback(self):
        """
        Submit feedback by showing a google form
        """
        webbrowser.open("https://goo.gl/forms/X9OB6AQSJSiKScBs2")

<<<<<<< HEAD
    def check_for_update(self):
        response = urllib2.urlopen('https://raw.githubusercontent.com/khaiquangnguyen/Chemics/master/APP_VERSION.html')
=======
<<<<<<< HEAD
    def updateSelection(self):
        response = urllib2.urlopen('http://khaiquangnguyen.github.io/chemics_update.html')
=======
    def check_for_update(self):
        response = urllib2.urlopen('https://raw.githubusercontent.com/khaiquangnguyen/Chemics/master/APP_VERSION.html')
>>>>>>> ce5b242c238e196dfe20cf2038e62c260c5c95bd
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5
        html = response.read()
        update = int(html[0])
        self.showUpdateDialog(update)

    def select_files(self):
        """
        Select the files which store the file
        """
        dialog = QFileDialog()
        files = dialog.getOpenFileNames()[0]
        if files:
            self.controller.reset()
            self.reset()
            self.controller.files = files
            self.controller.run()

<<<<<<< HEAD
    def move_progress_bar_forward(self, message=None, max_value=None, complete=False, value=1):
=======
<<<<<<< HEAD
    def makeProgress(self, message = None, maxValue = None, complete = False, value = 1):
=======
    def move_progress_bar_forward(self, message=None, max_value=None, complete=False, value=1):
>>>>>>> ce5b242c238e196dfe20cf2038e62c260c5c95bd
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5
        """
        Progress the progress dialog bar
        :param max_value: the maximum of the progress dialog. Used to signal reset progress dialog
        :param complete: whether the progress is the completion of the whole procedure
        :param message: the message to show on the progress dialog
        """

        if max_value is not None:
            self.progress = 0
            if self.progress_dialog is not None:
                self.progress_dialog.reset()
                self.progress_dialog.setValue(0)
                self.progress_dialog.setRange(0,max_value)
                if message is not None:
                    self.progress_dialog.setLabelText(message)
                self.progress_dialog.setWindowModality(Qt.WindowModal)
                self.progress_dialog.show()
            else:
<<<<<<< HEAD
=======
<<<<<<< HEAD
                self.progressDialog = QProgressDialog("Tasks in progress...", "Cancel", 0, maxValue, self)
                self.progressDialog.canceled.connect(self.cancelProgress)
                self.progressDialog.setWindowModality(Qt.WindowModal)
                self.progressDialog.show()

        else:
            if self.progressDialog is not None:
                if complete == True:
                    self.progressDialog.setValue(self.progressDialog.maximum())
                    self.progressDialog.reset()
=======
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5
                self.progress_dialog = QProgressDialog("Tasks in progress...", "Cancel", 0, max_value, self)
                self.progress_dialog.canceled.connect(self.cancel_progress_bar)
                self.progress_dialog.setWindowModality(Qt.WindowModal)
                self.progress_dialog.show()

        else:
            if self.progress_dialog is not None:
                if complete is True:
                    self.progress_dialog.setValue(self.progress_dialog.maximum())
                    self.progress_dialog.reset()
<<<<<<< HEAD
=======
>>>>>>> ce5b242c238e196dfe20cf2038e62c260c5c95bd
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5
                else:
                    if message is not None:
                        self.progress_dialog.setLabelText(message)
                    self.progress += 1
                    self.progress_dialog.setValue(self.progress)
                qApp.processEvents()

<<<<<<< HEAD
    def show_error_dialog(self, error_message='Unknown Error!'):
=======
<<<<<<< HEAD
    def showError(self, errorMessage = 'Unknown Error!'):
=======
    def show_error_dialog(self, error_message='Unknown Error!'):
>>>>>>> ce5b242c238e196dfe20cf2038e62c260c5c95bd
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5
        """
        Show the error message
        :param errorMessage: The message to show in the error message
        """
        if self.progress_dialog is not None:
            self.progress_dialog.reset()
            self.progress_dialog = None
        warning = QMessageBox()
        warning.setIcon(QMessageBox.Warning)
        warning.setText(error_message)
        warning.exec_()

<<<<<<< HEAD
=======
<<<<<<< HEAD
    def cancelProgress(self):
        """
        Action when the cancel button of the progress bar is clicked
        """
        self.controller.cancelProgress()

=======
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5
    def cancel_progress_bar(self):
        """
        Action when the cancelling_progress_bar button of the progress bar is clicked
        """
        self.controller.cancel_progress_bar()
<<<<<<< HEAD
=======
>>>>>>> ce5b242c238e196dfe20cf2038e62c260c5c95bd
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5

    def resizeEvent(self, resizeEvent):
        self.centralWidget().resize()

<<<<<<< HEAD
=======
<<<<<<< HEAD
    def update_dp_dnlog_figures(self, adjustedFigure, diaFigure):
        self.centralWidget().graphWidget.dpAndDnlogView.dpView.updateFigure(adjustedFigure)
        self.centralWidget().graphWidget.dpAndDnlogView.dNlogView.updateFigure(diaFigure)

    def updateTempOrMinFigure(self, aFigure):
        self.centralWidget().graphWidget.tempAndMinView.updateFigure(aFigure)

    def calKappa(self):
        # if self.controller.completedStep >=2:
        self.controller.calKappa()
        self.centralWidget().switchToKappa()
        self.controller.makeKappaGraph()

    def updateGeneralInfo(self):
        self.centralWidget().infoWidget.infoTable.updateGeneralInfo()

    def updateBasicPeakInfo(self):
        self.centralWidget().infoWidget.infoTable.updateBasicPeakInfo()

    def updateSigFitPeakInfo(self):
        self.centralWidget().infoWidget.infoTable.updateSigFitPeakInfo()

    def reset(self):
        self.centralWidget().switchToPeak()
=======
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5
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
            if (self.controller.is_usable_for_sigmoid_fit_list[curr_scan]):
                self.centralWidget().graph_widget.buttons_widget.change_scan_status_button.setText("Disable Scan")
            else:
                self.centralWidget().graph_widget.buttons_widget.change_scan_status_button.setText("Enable Scan")

    def update_scan_information_after_sigmoid_fit(self):
        self.centralWidget().info_widget.information_table.update_scan_information_after_sigmoid_fit()
        curr_scan = self.controller.current_scan
        if curr_scan in self.controller.unfinished_sigmoid_fit_scans_list:
            self.centralWidget().graph_widget.buttons_widget.change_scan_status_button.setText("Disable Scan")
            return
        if self.controller.min_pos_SMPS_list[curr_scan] is None or self.controller.min_pos_CCNC_list[curr_scan] is None:
            self.centralWidget().graph_widget.buttons_widget.change_scan_status_button.setText("Enable Scan")
        else:
            if self.controller.is_usable_for_kappa_cal_list[curr_scan] and self.controller.is_usable_for_sigmoid_fit_list[
                                                                                                            curr_scan]:
                self.centralWidget().graph_widget.buttons_widget.change_scan_status_button.setText("Disable Scan")
            else:
                self.centralWidget().graph_widget.buttons_widget.change_scan_status_button.setText("Enable Scan")

    def reset(self):
        self.centralWidget().switch_to_scan_widget()
<<<<<<< HEAD

    def update_kappa_values(self,sigma,temp,dd1,i1,dd2,i2,solu):
=======
>>>>>>> ce5b242c238e196dfe20cf2038e62c260c5c95bd

    def updateKappaVars(self,sigma,temp,dd1,i1,dd2,i2,solu):
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5
        self.controller.sigma = sigma
        self.controller.temp = temp
        self.controller.dd = dd1
        self.controller.iKappa = i1
        self.controller.dd2 = dd2
        self.controller.iKappa2 = i2
        self.controller.solubility = solu

<<<<<<< HEAD
=======
<<<<<<< HEAD
    def updateKappaGraph(self):
        self.centralWidget().graphWidget.graphView.updateFigure(self.controller.kappaGraph)

    def InputFlowRate(self):
        returnValue = self.controller.flowRate
        while True:
            input = QInputDialog.getDouble(self, self.tr("Get Flow Rate"),self.tr("Q(flow rate)"),0.3)
            if input[1] == True:
                returnValue = float(input[0])
                break
        return returnValue
=======
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5
    def update_kappa_info_and_graph(self):
        self.centralWidget().graph_widget.update_figure(self.controller.kappa_figure)
        self.centralWidget().info_widget.update_data()
        self.centralWidget().resize()
        self.centralWidget().info_widget.information_table.toggle_color(self.controller.current_kappa_point_index)
        ss = self.controller.kappa_points_data_list[self.controller.current_kappa_point_index][1]
        dp = self.controller.kappa_points_data_list[self.controller.current_kappa_point_index][0]
        if self.controller.is_show_all_k_points:
            status = self.controller.kappa_points_is_included_list[(dp, ss)]
            if status:
                self.centralWidget().graph_widget.control_widget.toggle_k_point_status_button.setText("Disable")
            else:
                self.centralWidget().graph_widget.control_widget.toggle_k_point_status_button.setText("Enable")
        self.centralWidget().setFocus()

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
<<<<<<< HEAD
=======
>>>>>>> ce5b242c238e196dfe20cf2038e62c260c5c95bd
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5

class ControlPanel(QWidget):
    def __init__(self, main_window=None):
        QWidget.__init__(self)
        # addSubWidget
        self.main_window = main_window
        self.layout = QHBoxLayout()
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
<<<<<<< HEAD
=======
<<<<<<< HEAD
        self.infoWidget = PeakTextDataWidget(self.mainWindow)
        self.graphWidget = PeakGraphWidget(self.mainWindow)
        self.layout.addWidget(self.infoWidget)
        self.layout.addWidget(self.graphWidget)
        self.setLayout(self.layout)

    def resize(self):
        self.infoWidget.resize(self.width(), self.height())
        self.graphWidget.resize(self.width(),self.height())

    def switchToKappa(self):
        self.clearLayout(self.layout)
        self.infoWidget = KappaTextDataWidget(self.mainWindow)
        self.graphWidget = KappaGraphWidget(self.mainWindow)
        self.clearLayout(self.layout)
        self.layout.addWidget(self.infoWidget)
        self.layout.addWidget(self.graphWidget)
        self.setLayout(self.layout)
        self.infoWidget.updateData()
        self.resize()

    def switchToPeak(self):
        self.clearLayout(self.layout)
        self.infoWidget = PeakTextDataWidget(self.mainWindow)
        self.graphWidget = PeakGraphWidget(self.mainWindow)
        self.layout.addWidget(self.infoWidget)
        self.layout.addWidget(self.graphWidget)
=======
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5
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
<<<<<<< HEAD
=======
>>>>>>> ce5b242c238e196dfe20cf2038e62c260c5c95bd
>>>>>>> 5ef947266cae82fc465cb122c098c082febbf7b5
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
        # finish kappa
        if self.main_window.controller.kappa_ax:
            self.main_window.controller.on_key_release_kappa_graph(event)
        else:
            self.main_window.controller.on_key_release(event)
