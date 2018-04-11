###############################
#
# IMPORT STATEMENTS
#
###############################


import socket
import sys
import urllib2
import Controller
import matplotlib
import webbrowser
from PySide.QtCore import *
from AlignmentUI import *
from KappaUI import *
from Graphs import *
# matplotlib.rcParams['backend.qt4'] = 'PySide'
###############################
#
# VARIABLES
#
###############################

class MainView(QMainWindow):
    """
    The main window (view) of the program
    """

    def __init__(self):
        QMainWindow.__init__(self)
        # create controller
        self.controller = Controller.Controller(self)
        self.progress_dialog = None
        self.setWindowTitle('Chemics')
        self.setMinimumHeight(800)
        self.setMinimumWidth(800)
        # set style
        # QApplication.setStyle(QStyleFactory.create('Cleanlooks'))
        self.raw_conc_time_graph = ConcOverTimeRawDataGraph()
        self.smoothed_conc_time_graph = ConcOverTimeSmoothGraph()
        self.temp_graph = TemperatureGraph()
        self.ratio_dp_graph = RatioOverDiameterGraph()

        print self.centralWidget()
        # self.setCentralWidget(ControlPanel(self))
        # self.setCentralWidget(SquareFigureCanvas())
        self.showMaximized()

        # create progress bar
        self.create_progress_bar()
        # create menu
        self.create_menus()
        # create central widgets
        self.central_widget_alignment = None
        self.create_central_widget()
        # create left dock widget for information related stuff
        self.create_left_docker()

    def create_progress_bar(self):
        self.progress_dialog = QProgressDialog("Starting..", "Cancel", 0, 100, self)
        # self.progress_dialog.canceled.connect(self.cancel_progress_bar)
        self.progress_dialog.setWindowModality(Qt.WindowModal)
        self.progress_dialog.setAutoReset(True)
        self.progress_dialog.setAutoClose(True)

    def create_central_widget(self):
        # create alignment central widget
        self.central_widget_alignment = CentralWidgetAlignment(self, self.raw_conc_time_graph,
                                                               self.smoothed_conc_time_graph, self.ratio_dp_graph,
                                                               self.temp_graph)
        self.setCentralWidget(self.central_widget_alignment)
        # self.central_widget_sigmoid_fit = CentralWidgetSigmoidFit()
        # self.central_widget_kappa = CentralWidgetKappaCalc()

    def create_menus(self):
        # Add navigation menu
        # Create the import file button
        select_action = QAction('&Select Files', self, shortcut="Ctrl+O", triggered =
        self.select_files)
        feedback_action = QAction('&Feedback', self, triggered=self.submit_feedback)
        self.menuBar().addAction(select_action)
        self.menuBar().addAction(feedback_action)

    def create_left_docker(self):
        left_widget = DockerAlignment("Information", self)
        left_widget.setAllowedAreas(Qt.RightDockWidgetArea|Qt.LeftDockWidgetArea)
        left_widget.setFeatures(QDockWidget.DockWidgetMovable | QDockWidget.DockWidgetClosable)
        self.menuBar().addAction(left_widget.toggleViewAction())
        self.addDockWidget(Qt.LeftDockWidgetArea, left_widget)

    def select_files(self):
        """
        Show a GUI to select the files which store the data.
        """
        files = QFileDialog.getOpenFileNames(self,'Open file', '', "Data files (*.csv *.txt)")[0]
        if files:
            # read in new files
            self.controller.start_or_restart(files)

    def get_counts_to_conc_conv(self):
        while True:
            input = QInputDialog.getDouble(self, "Get Counts2ConcConv", "Get Counts2ConcConv", 0.3, decimals=2)
            if input[1]:
                return float(input[0])

    @Slot(str)
    def init_progress_bar(self,task_name):
        self.progress_dialog.setLabelText(task_name)
        self.progress_dialog.setValue(0)
        self.progress_dialog.show()

    @Slot(int)
    def update_progress_bar(self,percentage):
        self.progress_dialog.setValue(percentage)
        qApp.processEvents()

    @Slot()
    def close_progress_bar(self):
        self.progress_dialog.reset()
        self.controller.cancel_progress_bar()

    def show_error_dialog(self, error_message='Unknown Error!'):
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

    def show_alignment_dialog(self):
        dialog = QMessageBox()
        dialog.setText("Select your preferred way to align the data")
        dialog.setInformativeText("Do you want to manually align the data, or let the program do it for you?")
        dialog.setIcon(QMessageBox.Question)
        manual_button = dialog.addButton("Manual",QMessageBox.RejectRole)
        auto_button = dialog.addButton("Auto",QMessageBox.AcceptRole)
        dialog.exec_()
        if dialog.clickedButton() == auto_button:
            next_dialog = SmoothAlgoDialog(self.controller.scans,self.controller)
            next_dialog.exec_()


    def calculate_all_kappa_values(self):
        """
        calculate the kappa value
        :return:
        """
        self.centralWidget().switch_to_kappa_widget()
        self.controller.calculate_all_kappa_values()
        self.controller.update_kappa_info_and_graph()

    def update_experiment_information(self):
        """
        update the general information of the entire dataset
        :return:
        """
        self.centralWidget().info_widget.information_table.update_experiment_information()

    def update_scan_information(self):
        """
        update the specific information of each scan at the beginning of the program
        :return:
        """
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
        """
        update the specific information of each scan after sigmoid fit
        :return:
        """
        self.centralWidget().info_widget.information_table.update_scan_information_after_sigmoid_fit()
        curr_scan = self.controller.current_scan
        if curr_scan in self.controller.unfinished_sigmoid_fit_scans_list:
            self.centralWidget().graph_widget.buttons_widget.change_scan_status_button.setText("Disable Scan")
            return
        if self.controller.min_pos_SMPS_list[curr_scan] is None or self.controller.min_pos_CCNC_list[curr_scan] is None:
            self.centralWidget().graph_widget.buttons_widget.change_scan_status_button.setText("Enable Scan")
        else:
            if self.controller.is_usable_for_kappa_cal_list[curr_scan] and \
                    self.controller.is_usable_for_sigmoid_fit_list[
                        curr_scan]:
                self.centralWidget().graph_widget.buttons_widget.change_scan_status_button.setText("Disable Scan")
            else:
                self.centralWidget().graph_widget.buttons_widget.change_scan_status_button.setText("Enable Scan")

    def reset(self):
        """
        reset the entire view
        :return:
        """
        # self.centralWidget().switch_to_scan_widget()

    def update_kappa_values(self, sigma, temp, dd1, i1, dd2, i2, solu):
        """
        update the kappa values
        :param sigma:
        :param temp:
        :param dd1:
        :param i1:
        :param dd2:
        :param i2:
        :param solu:
        :return:
        """
        self.controller.sigma = sigma
        self.controller.temp = temp
        self.controller.dd = dd1
        self.controller.iKappa = i1
        self.controller.dd2 = dd2
        self.controller.iKappa2 = i2
        self.controller.solubility = solu

    def update_kappa_info_and_graph(self):
        """
        update the kappa information and change the graph
        :return:
        """
        self.centralWidget().graph_widget.update_figure(self.controller.kappa_figure)
        self.centralWidget().info_widget.update()
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

    def submit_feedback(self):
        """
        Submit feedback by showing a google form
        """
        webbrowser.open("https://goo.gl/forms/X9OB6AQSJSiKScBs2")

    def update_alignment_graphs(self, a_scan):
        self.update_raw_conc_time_graph(a_scan)
        self.update_smooth_conc_time_graph(a_scan)
        self.update_temp_graph(a_scan)
        self.update_ratio_dp_graph(a_scan)

    def update_raw_conc_time_graph(self, a_scan):
        self.raw_conc_time_graph.update_graph(a_scan)

    def update_smooth_conc_time_graph(self, a_scan):
        self.smoothed_conc_time_graph.update_graph(a_scan)

    def update_temp_graph(self, a_scan):
        self.temp_graph.update_graph(a_scan)

    def update_ratio_dp_graph(self, a_scan):
        self.ratio_dp_graph.update_graph(a_scan)

#######################################################
#
# CONTROL PANEL CLASS. CONTAINS ALL THE CONTROL BUTTONS OF THE UI
#
#######################################################

class CentralWidgetAlignment(QWidget):
    def __init__(self,parent,raw_conc_time_graph, smoothed_conc_time_graph, ratio_dp_graph, temp_graph):
        super(self.__class__, self).__init__()

        # Add widgets
        # init the necessary contents
        self.raw_conc_time_graph = raw_conc_time_graph
        self.smoothed_conc_time_graph = smoothed_conc_time_graph
        self.ratio_dp_graph = ratio_dp_graph
        self.temp_graph = temp_graph
        self.h_splitter_1 = QSplitter(Qt.Horizontal)
        self.h_splitter_1.addWidget(self.raw_conc_time_graph)
        self.h_splitter_1.addWidget(self.smoothed_conc_time_graph)
        self.h_splitter_2 = QSplitter(Qt.Horizontal)
        self.h_splitter_2.addWidget(self.ratio_dp_graph)
        self.h_splitter_2.addWidget(self.temp_graph)
        self.v_splitter = QSplitter(Qt.Vertical)
        self.v_splitter.addWidget(self.h_splitter_1)
        self.v_splitter.addWidget(self.h_splitter_2)
        hbox = QHBoxLayout(self)
        hbox.addWidget(self.v_splitter)
        self.setLayout(hbox)


class DockerAlignment(QDockWidget):
    def __init__(self,title, parent=None):
        super(self.__class__, self).__init__(title,parent)
        form_layout = QFormLayout()
        h_layout = QHBoxLayout()
        h_layout.addWidget()




if __name__ == "__main__":
    app = QApplication(sys.argv)
    main_window = MainView()
    main_window.show()
    sys.exit(app.exec_())
