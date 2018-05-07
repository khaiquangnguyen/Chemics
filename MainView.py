import webbrowser

from Graphs import *
from HelperFunctions import *
from CustomCentralWidgets import *
from CustomDockerWidgets import *
from CustomModalDialogs import *
from Controller import Controller
import gc

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
        self.controller = Controller(self)
        self.setWindowTitle('Chemics')
        # set style
        self.raw_conc_time_graph = ConcOverTimeRawDataGraph()
        self.smoothed_conc_time_graph = ConcOverTimeSmoothGraph()
        self.temp_graph = TemperatureGraph()
        self.ratio_dp_graph = RatioOverDiameterGraph()
        self.kappa_graph = KappaGraph()
        # the window menu.
        # create central widgets
        self.showMaximized()
        # create progress bar
        self.progress_dialog = None
        self.create_progress_bar()
        # create menu
        self.create_menus()
        # create central widget
        self.stacked_central_widget = None
        self.central_widget_alignment = None
        self.central_widget_kappa = None
        self.create_central_widget()
        # create left dock widget for information related stuff
        self.align_docker_widget = None
        self.create_align_docker()
        self.setDockOptions(QMainWindow.VerticalTabs | QMainWindow.AnimatedDocks | QMainWindow.ForceTabbedDocks)
        # init the right menu
        self.set_menu_bar_by_stage()
        # some other cosmetics
        self.font = QFont("Calibri")
        app.setFont(self.font)

    def create_progress_bar(self):
        self.progress_dialog = QProgressDialog("Starting..", "Cancel", 0, 100, self)
        self.progress_dialog.setWindowModality(Qt.WindowModal)
        self.progress_dialog.setAutoReset(True)
        self.progress_dialog.setAutoClose(True)

    def create_central_widget(self):
        # create alignment central widget
        self.central_widget_alignment = CentralWidgetAlignment(self, self.raw_conc_time_graph,
                                                               self.smoothed_conc_time_graph, self.ratio_dp_graph,
                                                               self.temp_graph)
        self.central_widget_kappa = CentralWidgetKappa(self, self.kappa_graph)
        self.stacked_central_widget = QStackedWidget()
        self.stacked_central_widget.addWidget(self.central_widget_alignment)
        self.stacked_central_widget.addWidget(self.central_widget_kappa)
        # lock out the menus that we will not use
        self.setCentralWidget(self.stacked_central_widget)
        # Lock down all actions in the action menu
        self.action_menu.setDisabled(True)

    def create_menus(self):
        # Add file menu
        self.file_menu = QMenu("&File")
        new_action = QAction('&New Project from Files...', self, shortcut="Ctrl+N", triggered=self.open_files)
        open_action = QAction('&Open Existing Project...', self, shortcut="Ctrl+O", triggered=self.open_project)
        save_action = QAction('&Save Project', self, shortcut="Ctrl+S", triggered=self.save_project)
        save_as_action = QAction('&Save Project As', self, triggered=self.save_project_as)
        reset_action = QAction('&Reset Project', self, triggered=self.reset_project)
        export_data_action = QAction('&Export Project Data', self, triggered=self.reset_project)
        exit_action = QAction('&Exit', self, triggered=app.quit)
        self.file_menu.addAction(new_action)
        self.file_menu.addSeparator()
        self.file_menu.addActions([open_action, save_action, save_as_action, export_data_action])
        self.file_menu.addSeparator()
        self.file_menu.addAction(reset_action)
        self.file_menu.addSeparator()
        self.file_menu.addAction(exit_action)
        self.menuBar().addMenu(self.file_menu)
        # add action menu
        self.action_menu = QMenu("&Actions")
        preview_all_action = QAction('Preview all scans', self, triggered=self.preview_all_scans)
        smooth_data_action = QAction('&Smooth Data', self, triggered=self.show_smooth_dialog)
        auto_align_action = QAction('&Auto Align', self, triggered=self.show_auto_align_dialog)
        correct_charges = QAction('Correct Charges All Scans', self, triggered=self.correct_charges)
        correct_charges_one = QAction('Correct Charges One Scan', self, triggered=self.correct_charges_one)
        auto_fit_action = QAction('&Auto Fit Sigmoid', self, triggered=self.show_auto_fit_sigmoid_dialog)
        cal_kappa_action = QAction('&Calculate Kappa', self, triggered=self.show_kappa_params_dialog)
        self.action_menu.addAction(preview_all_action)
        self.action_menu.addSeparator()
        self.action_menu.addActions([smooth_data_action, auto_align_action])
        self.action_menu.addSeparator()
        self.action_menu.addActions([correct_charges, correct_charges_one, auto_fit_action])
        self.action_menu.addSeparator()
        self.action_menu.addAction(cal_kappa_action)
        self.menuBar().addMenu(self.action_menu)
        # add settings action
        setting_action = QAction('&Settings', self, triggered=self.show_setting_dialog)
        self.menuBar().addAction(setting_action)
        # add window menu
        self.window_menu = QMenu("&Windows")
        self.menuBar().addMenu(self.window_menu)
        # add Help menu
        self.help_menu = QMenu("&Help")
        feedback_action = QAction('&Send Feedback', self, triggered=self.submit_feedback)
        check_for_update_action = QAction('&Check for Updates', self, triggered=self.submit_feedback)
        contact_creator_action = QAction('&Contact Creator', self, triggered=self.submit_feedback)
        user_manual_action = QAction('&User Manual', self, triggered=self.submit_feedback)
        self.help_menu.addActions(
            [feedback_action, check_for_update_action, contact_creator_action, user_manual_action])
        self.menuBar().addMenu(self.help_menu)

    def set_menu_bar_by_stage(self):
        if self.controller.stage == "init":
            self.action_menu.setDisabled(True)
            self.window_menu.actions()[1].setDisabled(True)
            self.window_menu.actions()[2].setDisabled(True)
            self.window_menu.actions()[3].setDisabled(True)
            file_action_list = self.file_menu.actions()
            file_action_list[3].setDisabled(True)
            file_action_list[4].setDisabled(True)
            file_action_list[5].setDisabled(True)
            file_action_list[6].setDisabled(True)
            file_action_list[7].setDisabled(True)
        elif self.controller.stage == "align":
            # Enable action menu and lock down the ability to edit sigmoid
            self.window_menu.actions()[1].setDisabled(True)
            self.window_menu.actions()[2].setDisabled(True)
            self.window_menu.actions()[3].setDisabled(True)
            self.action_menu.setDisabled(False)
            action_list = self.action_menu.actions()
            action_list[6].setDisabled(True)
            action_list[7].setDisabled(True)
            action_list[8].setDisabled(True)
            action_list[9].setDisabled(True)
            file_action_list = self.file_menu.actions()
            file_action_list[3].setDisabled(False)
            file_action_list[4].setDisabled(False)
            file_action_list[5].setDisabled(False)
            file_action_list[6].setDisabled(False)
            file_action_list[7].setDisabled(False)
        elif self.controller.stage == "sigmoid":
            # Edit action menu
            action_list = self.action_menu.actions()
            self.action_menu.setDisabled(False)
            self.window_menu.actions()[1].setDisabled(False)
            self.window_menu.actions()[2].setDisabled(True)
            self.window_menu.actions()[3].setDisabled(True)
            action_list[6].setEnabled(True)
            action_list[7].setEnabled(True)
            action_list[8].setEnabled(True)
            action_list[9].setEnabled(True)
            file_action_list = self.file_menu.actions()
            # enable all file action
            for i in range(len(file_action_list)):
                file_action_list[i].setDisabled(False)
        elif self.controller.stage == "kappa":
            self.window_menu.actions()[1].setDisabled(False)
            self.window_menu.actions()[2].setDisabled(False)
            self.window_menu.actions()[3].setDisabled(False)

            self.action_menu.setDisabled(False)
            action_list = self.action_menu.actions()
            action_list[6].setDisabled(True)
            action_list[7].setDisabled(True)
            action_list[8].setDisabled(True)
            action_list[9].setDisabled(True)
            file_action_list = self.file_menu.actions()
            file_action_list[3].setDisabled(False)
            file_action_list[4].setDisabled(False)
            file_action_list[5].setDisabled(False)
            file_action_list[6].setDisabled(False)
            file_action_list[7].setDisabled(False)
            action_list[6].setEnabled(True)
            action_list[7].setEnabled(True)
            action_list[8].setEnabled(True)
            action_list[9].setEnabled(True)
            self.window_menu.actions()[1].setDisabled(False)
            self.window_menu.actions()[2].setDisabled(False)

    def show_setting_dialog(self):
        setting_dialog = SettingDialog(self)
        setting_dialog.exec_()

    def create_align_docker(self):
        # create alignment docker
        self.align_docker = QDockWidget("Scan Information", self)
        self.align_docker_widget = DockerWidgetAlignment(self.controller)
        self.align_docker.setWidget(self.align_docker_widget)
        self.align_docker.setAllowedAreas(Qt.RightDockWidgetArea | Qt.LeftDockWidgetArea)
        self.align_docker.setFeatures(QDockWidget.DockWidgetMovable | QDockWidget.DockWidgetClosable)
        self.window_menu.addAction(self.align_docker.toggleViewAction())
        self.addDockWidget(Qt.LeftDockWidgetArea, self.align_docker)
        # create sigmoid docker
        self.sigmoid_docker = QDockWidget("Sigmoid Parameters", self)
        self.sigmoid_docker_widget = DockerSigmoidWidget(self.controller)
        self.sigmoid_docker.setWidget(self.sigmoid_docker_widget)
        self.sigmoid_docker.setAllowedAreas(Qt.RightDockWidgetArea | Qt.LeftDockWidgetArea)
        self.sigmoid_docker.setFeatures(QDockWidget.DockWidgetMovable | QDockWidget.DockWidgetClosable)
        self.window_menu.addAction(self.sigmoid_docker.toggleViewAction())
        self.addDockWidget(Qt.LeftDockWidgetArea, self.sigmoid_docker)
        self.tabifyDockWidget(self.sigmoid_docker, self.align_docker)
        # create kappa docker
        self.kappa_docker = QDockWidget("Kappa Values", self)
        self.kappa_docker_widget = DockerKappaWidget(self.controller, self.kappa_graph)
        self.kappa_docker.setWidget(self.kappa_docker_widget)
        self.kappa_docker.setAllowedAreas(Qt.RightDockWidgetArea | Qt.LeftDockWidgetArea)
        self.kappa_docker.setFeatures(QDockWidget.DockWidgetMovable | QDockWidget.DockWidgetClosable)
        self.window_menu.addAction(self.kappa_docker.toggleViewAction())
        self.addDockWidget(Qt.LeftDockWidgetArea, self.kappa_docker)
        show_window_action = QAction("&Show Alignment Graphs", self, triggered=self.switch_central_widget)
        show_window_action.setCheckable(True)
        self.window_menu.addAction(show_window_action)

    def show_sigmoid_docker(self):
        # self.addDockWidget(Qt.LeftDockWidgetArea, self.sigmoid_docker)
        # self.sigmoid_docker.show()
        # self.sigmoid_docker.raise_()
        # first, let's make them both unchecked
        if self.window_menu.actions()[0].isChecked():
            self.window_menu.actions()[0].trigger()
        if self.window_menu.actions()[1].isChecked():
            self.window_menu.actions()[1].trigger()
        if self.window_menu.actions()[2].isChecked():
            self.window_menu.actions()[2].trigger()
        # now, enable them one by one
        self.window_menu.actions()[1].trigger()
        self.window_menu.actions()[0].trigger()

    def open_files(self):
        files = QFileDialog.getOpenFileNames(self, "Open file", '', "Data files (*.csv *.txt)")[0]
        if files:
            # read in new files
            self.controller.start(files)

    def open_project(self):
        project_file = QFileDialog.getOpenFileName(self, "Open file", '', "Project files (*.chemics)")[0]
        if project_file:
            # read in new files
            self.controller.load_project(project_file)

    def save_project(self):
        self.controller.save_project()

    def reset_project(self):
        self.controller.reset_project()

    def export_project_data(self):
        self.controller.export_project_data()

    def save_project_as(self):
        project_file = QFileDialog.getSaveFileName(self, "Save file", '', "Project files (*.chemics)")[0]
        if project_file:
            # read in new files
            self.controller.set_save_name(project_file)
            self.controller.save_project()

    def get_counts_to_conc_conv(self):
        cc = QInputDialog.getDouble(self, "Counts2ConcConv", "Counts2ConcConv", decimals=2)
        if cc[1]:
            return float(cc[0])
        else:
            return None

    def show_auto_or_manual_question_dialog(self):
        dialog = QMessageBox()
        dialog.setText("Select your preferred way to align the data")
        dialog.setInformativeText("Do you want to manually align the data, or let the program do it for you?")
        dialog.setIcon(QMessageBox.Question)
        manual_button = dialog.addButton("Manual", QMessageBox.RejectRole)
        auto_button = dialog.addButton("Auto", QMessageBox.AcceptRole)
        dialog.exec_()
        if dialog.clickedButton() == auto_button:
            smooth_dialog = SmoothAlgoDialog(self.controller)
            if smooth_dialog.exec_() == 1:
                next_dialog = SetBaseShiftDialog(self.controller)
                next_dialog.exec_()
        elif dialog.clickedButton() == manual_button:
            self.controller.prepare_data_for_manual_inputs()

    def preview_all_scans(self):
        timer, ok = QInputDialog.getDouble(self, "Time of each preview", "Enter the amount of pause time between each "
                                                                         "scan", value=0.5,
                                                                         minValue=0.1, maxValue=3, decimals=1)
        if timer and ok:
            self.controller.preview_scans(timer)

    def correct_charges(self):
        if len(self.controller.scans) == 0:
            self.show_error_by_type("no_data")
            return
        message = QMessageBox()
        message.setWindowTitle("Correct Charges All Scans")
        message.setIcon(QMessageBox.Question)
        message.setText("Are you sure you want to start correcting charges for all scans?")
        message.setInformativeText(
            "It is highly recommended that all scans are aligned before performing charge corrections!")
        message.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        message.setDefaultButton(QMessageBox.No)
        ret = message.exec_()
        if ret == QMessageBox.Yes:
            self.controller.correct_charges()
            self.show_auto_manual_fit_sigmoid_dialog()

    def correct_charges_one(self):
        if len(self.controller.scans) == 0:
            self.show_error_by_type("no_data")
            return
        message = QMessageBox()
        message.setWindowTitle("Correct Charges One Scan")
        message.setIcon(QMessageBox.Question)
        message.setText("Please confirm that you wish to correct charges for this scan!")
        message.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        message.setDefaultButton(QMessageBox.No)
        ret = message.exec_()
        if ret == QMessageBox.Yes:
            curr_scan = self.controller.scans[self.controller.curr_scan_index]
            curr_scan.generate_processed_data()
            curr_scan.correct_charges()
            self.controller.switch_to_scan(self.controller.curr_scan_index)

    def show_smooth_dialog(self):
        if len(self.controller.scans) == 0:
            self.show_error_by_type("no_data")
            return
        smooth_dialog = SmoothAlgoDialog(self.controller)
        smooth_dialog.exec_()

    def show_auto_align_dialog(self):
        if len(self.controller.scans) == 0:
            self.show_error_by_type("no_data")
            return
        align_dlg = SetBaseShiftDialog(self.controller)
        align_dlg.exec_()

    def show_auto_manual_fit_sigmoid_dialog(self):
        dialog = QMessageBox()
        dialog.setText("Select your preferred way to fit the sigmoid line to the data")
        dialog.setInformativeText("Do you want to manually fit the sigmoid lines, or let the program do it for you?")
        dialog.setIcon(QMessageBox.Question)
        manual_button = dialog.addButton("Manual", QMessageBox.RejectRole)
        auto_button = dialog.addButton("Auto", QMessageBox.AcceptRole)
        dialog.exec_()
        if dialog.clickedButton() == auto_button:
            self.controller.auto_fit_sigmoid()

    def show_auto_fit_sigmoid_dialog(self):
        dialog = QMessageBox()
        dialog.setText("Are you sure you want to let the program fit the sigmoid lines for you?")
        dialog.setIcon(QMessageBox.Question)
        dialog.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        dialog.setDefaultButton(QMessageBox.No)
        ret = dialog.exec_()
        if ret == QMessageBox.Yes:
            self.controller.auto_fit_sigmoid()

    def update_scan_info_and_graphs(self):
        a_scan = self.controller.scans[self.controller.curr_scan_index]
        self.raw_conc_time_graph.update_graph(a_scan)
        self.smoothed_conc_time_graph.update_graph(a_scan)
        self.temp_graph.update_graph(a_scan)
        self.ratio_dp_graph.update_graph(a_scan)
        self.align_docker_widget.update_scan_info()
        self.sigmoid_docker_widget.update_scan_info()

    def update_experiment_info(self):
        self.align_docker_widget.update_experiment_info()

    def show_kappa_params_dialog(self):
        kappa_dialog = SelectParamsKappaDialog(self.controller)
        kappa_dialog.exec_()

    def reset(self):
        """
        reset the entire view
        :return:
        """
        # self.centralWidget().switch_to_scan_widget()

    def switch_to_kappa_view(self):
        if self.window_menu.actions()[0].isChecked():
            self.window_menu.actions()[0].trigger()
        if self.window_menu.actions()[1].isChecked():
            self.window_menu.actions()[1].trigger()
        self.stacked_central_widget.setCurrentWidget(self.central_widget_kappa)
        self.kappa_graph.update_all_kappa_points(self.controller.alpha_pinene_dict,
                                                 self.controller.is_valid_kappa_points)
        self.kappa_docker_widget.update_kappa_values()
        self.set_menu_bar_by_stage()


        if not self.window_menu.actions()[2].isChecked():
            self.window_menu.actions()[2].trigger()

    def update_kappa_graph(self):
        self.kappa_docker_widget.update_kappa_graph()

    def submit_feedback(self):
        """
        Submit feedback by showing a google form
        """
        webbrowser.open("https://goo.gl/forms/X9OB6AQSJSiKScBs2")

    @Slot(str)
    def init_progress_bar(self, task_name):
        self.progress_dialog.setLabelText(task_name)
        self.progress_dialog.setValue(0)
        self.progress_dialog.show()

    @Slot(int)
    def update_progress_bar(self, percentage):
        self.progress_dialog.setValue(percentage)
        # qApp.processEvents()

    @Slot()
    def close_progress_bar(self):
        self.progress_dialog.reset()

    @staticmethod
    def show_error_message(title=None, text=None, subtext=None):
        warning = QMessageBox()
        warning.setIcon(QMessageBox.Warning)
        warning.setWindowTitle(title)
        warning.setText(text)
        warning.setInformativeText(subtext)
        warning.exec_()

    def show_error_by_type(self, type):
        if type == "no_data":
            title = "Error!"
            text = "There is no scan data to perform this action!"
            subtext = "Please import scan data through File/New or import a project through File/Open first!"
            self.show_error_message(title, text, subtext)

    def set_font(self,font,size):
        self.font = QFont(font.family(),size)
        app.setFont(self.font)

    def reset(self):
        # enable align docker, and disable the other dockers
        if not self.window_menu.actions()[0].isChecked():
            self.window_menu.actions()[0].trigger()
        if self.window_menu.actions()[1].isChecked():
            self.window_menu.actions()[1].trigger()
        if self.window_menu.actions()[2].isChecked():
            self.window_menu.actions()[2].trigger()
        self.stacked_central_widget.setCurrentWidget(self.central_widget_alignment)

    def switch_central_widget(self):
        self.stacked_central_widget.setCurrentIndex(self.stacked_central_widget.count() - 1 -
                                                    self.stacked_central_widget.currentIndex())





if __name__ == "__main__":
    app = QApplication(sys.argv)
    main_window = MainView()
    main_window.show()
    sys.exit(app.exec_())
