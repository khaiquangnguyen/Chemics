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
from HelperFunctions import *
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
        docker = QDockWidget("Information and Settings", self)
        content_widget = DockerWidgetAlignment()
        docker.setWidget(content_widget)
        docker.setAllowedAreas(Qt.RightDockWidgetArea|Qt.LeftDockWidgetArea)
        docker.setFeatures(QDockWidget.DockWidgetMovable | QDockWidget.DockWidgetClosable)
        self.menuBar().addAction(docker.toggleViewAction())
        self.addDockWidget(Qt.LeftDockWidgetArea, docker)

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

class LongLabel(QLabel):
    def paintEvent( self, event ):
        painter = QPainter(self)
        metrics = QFontMetrics(self.font())
        elided  = metrics.elidedText(self.text(), Qt.ElideRight, self.width())
        painter.drawText(self.rect(), self.alignment(), elided)

class DockerWidgetAlignment(QWidget):
    def __init__(self):
        super(self.__class__, self).__init__()
        # set up the layout
        form_layout = QFormLayout()
        # The basic information of a scan
        # add a title
        self.directory = LongLabel("C://asdfasfasdfasghgfhghfghfhgfhgfhgfdfasdfasfdss")
        self.directory.setMaximumWidth(150)
        self.directory.setToolTip(self.directory.text())
        self.directory.setAlignment(Qt.AlignRight)
        form_layout.addRow("Directory", self.directory)
        # add a title
        form_layout.addRow(TitleHLine("Experiment Information"))
        # add the items
        # -- add date
        self.experiment_date = QLabel("10/6/2017")
        self.experiment_date.setAlignment(Qt.AlignRight)
        form_layout.addRow("Date (m/d/y)", self.experiment_date)
        # -- add smps duration
        self.scan_duration = QLabel("135")
        self.scan_duration.setAlignment(Qt.AlignRight)
        form_layout.addRow("Scan duration(s)", self.scan_duration)
        #-- add Counts2ConcConv
        self.counts_2_conc_conv = QLabel("0.06")
        self.counts_2_conc_conv.setAlignment(Qt.AlignRight)
        form_layout.addRow("Counts2ConcConv",self.counts_2_conc_conv)

        #-------------------
        # add a title
        form_layout.addRow(TitleHLine("Scan Information"))
        #-- add the scan selector
        self.scan_selector = CustomSpinBox()
        form_layout.addRow("Scan Index", self.scan_selector)
        #-- add the shift selector
        self.shift_selector = CustomSpinBox()
        form_layout.addRow("Shift (s)", self.shift_selector)
        # -- add the super saturation indicator
        self.super_saturation = QLabel("0.2")
        self.super_saturation.setAlignment(Qt.AlignRight)
        form_layout.addRow("Super Saturation (%)", self.super_saturation)

        # -- add the start time
        self.scan_start_time = QLabel("10:26:25")
        self.scan_start_time.setAlignment(Qt.AlignRight)
        form_layout.addRow("Start Time (h:m:s)", self.scan_start_time)
        # -- ad the end time
        self.scan_end_time = QLabel("10:24:25")
        self.scan_end_time.setAlignment(Qt.AlignRight)
        form_layout.addRow("End Time (h:m:s)", self.scan_end_time)
        # add the layout
        self.setLayout(form_layout)

class TitleHLine(QWidget):
    def __init__(self,title):
        super(self.__class__, self).__init__()
        h_layout = QHBoxLayout()
        h_layout.addWidget(QHLine())
        h_layout.addWidget(QLabel(title))
        h_layout.addWidget(QHLine())
        self.setLayout(h_layout)


class QHLine(QFrame):
    def __init__(self):
        super(QHLine, self).__init__()
        self.setFrameShape(QFrame.HLine)
        self.setFrameShadow(QFrame.Raised)
        self.setLineWidth(2)


class QVLine(QFrame):
    def __init__(self):
        super(QVLine, self).__init__()
        self.setFrameShape(QFrame.VLine)
        self.setFrameShadow(QFrame.Sunken)

class CustomSpinBox(QWidget):
    def __init__(self):
        super(self.__class__, self).__init__()
        h_layout = QHBoxLayout()
        prev_button = QToolButton()
        prev_button.setArrowType(Qt.LeftArrow)
        next_button = QToolButton()
        next_button.setArrowType(Qt.RightArrow)
        content_box = QSpinBox()
        content_box.setButtonSymbols(QAbstractSpinBox.NoButtons)
        content_box.setMinimumWidth(80)
        content_box.setAlignment(Qt.AlignCenter)
        h_layout.addWidget(prev_button)
        h_layout.addWidget(content_box)
        h_layout.addWidget(next_button)
        h_layout.setAlignment(Qt.AlignRight)
        h_layout.setContentsMargins(0,0,0,0)
        self.setLayout(h_layout)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyle("cleanlooks")
    font = QFont("Calibri")
    app.setFont(font)
    main_window = MainView()
    main_window.show()
    sys.exit(app.exec_())
