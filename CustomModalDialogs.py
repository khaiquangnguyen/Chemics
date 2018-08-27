from Graphs import *
from CustomWidgets import *
import settings
from copy import deepcopy
from datetime import *
from PyQt5.QtCore import *
from scipy import *


class ScanInformationDialog(QDialog):
    def __init__(self, scan):
        super(self.__class__, self).__init__()
        # initiate the layout
        form_layout = QFormLayout()
        # scan status
        self.additional_information = QTextEdit("Welcome to Chemics!")
        self.additional_information.setReadOnly(True)
        self.additional_information.setMaximumHeight(100)
        if scan.is_valid():
            self.scan_status = QLabel("VALID")
            self.scan_status.setStyleSheet("QWidget { background-color:None}")
            self.additional_information.setText("The scan shows no problem.")
        else:
            self.scan_status = QLabel("INVALID")
            self.scan_status.setStyleSheet("QWidget { color: white; background-color:red}")
            self.additional_information.setText(scan.decode_status_code())
        form_layout.addRow("Scan status", self.scan_status)
        form_layout.addRow("Status Info", self.additional_information)
        # times: start time, end time, duration
        start_time = date.strftime(scan.start_time, "%H:%M:%S")
        end_time = date.strftime(scan.end_time, "%H:%M:%S")
        duration = str(scan.duration)
        time_group_box = QGroupBox()
        h_layout = QHBoxLayout()
        start_time_box = QLineEdit(start_time)
        start_time_box.setReadOnly(True)
        end_time_box = QLineEdit(end_time)
        end_time_box.setReadOnly(True)
        duration_box = QLineEdit(duration)
        duration_box.setReadOnly(True)
        h_layout.addWidget(QLabel("Start time"))
        h_layout.addWidget(start_time_box)
        h_layout.addWidget(QLabel("End Time"))
        h_layout.addWidget(end_time_box)
        h_layout.addWidget(QLabel("Duration"))
        h_layout.addWidget(duration_box)
        time_group_box.setLayout(h_layout)
        form_layout.addRow("Times and duration", time_group_box)
        self.setLayout(form_layout)


class SettingDialog(QDialog):
    def __init__(self, main_view):
        super(self.__class__, self).__init__()
        self.main_view = main_view
        form_layout = QFormLayout()
        # Font stuffs
        h_layout = QHBoxLayout()
        self.font_selector = QFontComboBox()
        self.font_selector.setCurrentFont(QFont("Calibri"))
        self.font_selector.currentFontChanged.connect(self.fontChanged)
        h_layout.addWidget(self.font_selector)
        self.size_selector = QSpinBox()
        self.size_selector.setValue(12)
        self.size_selector.valueChanged.connect(self.fontSizeChanged)
        h_layout.addWidget(self.size_selector)
        form_layout.addRow("Fonts", h_layout)
        # # Matplotlib style
        # self.plot_style_combobox = QComboBox()
        # style_list = ['default', 'classic'] + sorted(
        #     style for style in plt.style.available if style != 'classic')
        # self.plot_style_combobox.addItems(style_list)
        # self.plot_style_combobox.currentIndexChanged.connect(self.plotStyleChanged)
        # form_layout.addRow("Plot Style", self.plot_style_combobox)
        button_boxes = QDialogButtonBox()
        apply_button = button_boxes.addButton("Ok", QDialogButtonBox.AcceptRole)
        apply_button.clicked.connect(self.accept)
        form_layout.addWidget(button_boxes)
        self.setLayout(form_layout)

    def fontChanged(self, font):
        self.main_view.set_font(font, self.size_selector.value())

    def fontSizeChanged(self, size):
        self.main_view.set_font(self.font_selector.currentFont(), size)

    def plotStyleChanged(self):
        style = self.plot_style_combobox.currentText()
        self.main_view.change_pltstyle(style)


class SelectParamsKappaDialog(QDialog):
    def __init__(self, controller):
        super(self.__class__, self).__init__()
        self.controller = controller
        self.setWindowTitle("Select parameters for kappa calculation!")
        self.sigma_spinbox = QDoubleSpinBox()
        self.sigma_spinbox.setValue(self.controller.sigma)
        self.sigma_spinbox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.temp_spinbox = QDoubleSpinBox()
        self.temp_spinbox.setValue(self.controller.temp)
        self.temp_spinbox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.dd_1_spinbox = QDoubleSpinBox()
        self.dd_1_spinbox.setValue(self.controller.dd_1)
        self.dd_1_spinbox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.dd_2_spinbox = QDoubleSpinBox()
        self.dd_2_spinbox.setValue(self.controller.dd_2)
        self.dd_2_spinbox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.i_kappa_1_spinbox = QDoubleSpinBox()
        self.i_kappa_1_spinbox.setValue(self.controller.i_kappa_1)
        self.i_kappa_1_spinbox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.i_kappa_2_spinbox = QDoubleSpinBox()
        self.i_kappa_2_spinbox.setValue(self.controller.i_kappa_2)
        self.i_kappa_2_spinbox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.solubility_spinbox = QDoubleSpinBox()
        self.solubility_spinbox.setValue(self.controller.solubility)
        self.solubility_spinbox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        form_layout = QFormLayout()
        form_layout.addRow(self.tr("&Sigma"), self.sigma_spinbox)
        form_layout.addRow(self.tr("&Temperature"), self.temp_spinbox)
        form_layout.addRow(self.tr("&dry diameter(1)"), self.dd_1_spinbox)
        form_layout.addRow(self.tr("&iKappa(1)"), self.i_kappa_1_spinbox)
        form_layout.addRow(self.tr("&dry diameter(2)"), self.dd_2_spinbox)
        form_layout.addRow(self.tr("&iKappa(2)"), self.i_kappa_2_spinbox)
        form_layout.addRow(self.tr("&solubility"), self.solubility_spinbox)
        button_boxes = QDialogButtonBox()
        apply_button = button_boxes.addButton("Apply", QDialogButtonBox.ApplyRole)
        apply_button.clicked.connect(self.apply)
        form_layout.addWidget(button_boxes)
        self.setLayout(form_layout)

    def apply(self):
        # transfer all of the values to the controller
        self.controller.set_sigma(self.sigma_spinbox.value())
        self.controller.set_temp(self.temp_spinbox.value())
        self.controller.set_dd_1(self.dd_1_spinbox.value())
        self.controller.set_dd_2(self.dd_2_spinbox.value())
        self.controller.set_i_kappa_1(self.i_kappa_1_spinbox.value())
        self.controller.set_i_kapp_2(self.i_kappa_2_spinbox.value())
        self.controller.set_solubility(self.solubility_spinbox.value())
        self.accept()
        self.controller.cal_kappa()


class SelectParamsSigmoidDialog(QDialog):
    def __init__(self, controller):
        super(self.__class__, self).__init__()
        # the master layout
        form_layout = QFormLayout()
        self.controller = controller
        self.num_sigmoid_lines = 0
        self.dp_widgets = []
        # set up the groupbox which contains the graph of the scan
        # we assume that a scan in the middle will look very nice.
        self.curr_scan_index = self.controller.curr_scan_index
        self.setWindowTitle("Select parameters for sigmoid fit")
        self.graph = RatioOverDiameterGraph(self)
        self.graph.update_graph(self.scans[self.curr_scan_index])
        form_layout.addRow(self.graph)
        # -- add the area to select params for sigmoid
        sigmoid_line_spinbox = ArrowSpinBox(forward=True)
        form_layout.addRow("Number of sigmoid lines", sigmoid_line_spinbox)
        sigmoid_line_spinbox.setCallback(self.num_params_changed)
        # -- add the space to contain the param inputs for sigmoid lines
        self.params_container = QVBoxLayout()
        # the final accept button
        form_layout.addRow(self.params_container)
        # -- add the final apply button
        button_boxes = QDialogButtonBox()
        apply_button = button_boxes.addButton("Apply", QDialogButtonBox.ApplyRole)
        apply_button.clicked.connect(self.apply)
        form_layout.addWidget(button_boxes)
        self.setLayout(form_layout)

    def num_params_changed(self, n):
        if n > self.num_sigmoid_lines:
            for i in range(n - self.num_sigmoid_lines):
                self.add_params_group_box()
        else:
            for i in range(self.num_sigmoid_lines - n):
                self.sub_params_group_box()

    def sub_params_group_box(self):
        self.num_sigmoid_lines = max(self.num_sigmoid_lines - 1, 0)
        if len(self.dp_widgets) > 0:
            to_del = self.params_container.takeAt(self.params_container.count() - 1)
            del self.dp_widgets[-1]
            to_del.widget().deleteLater()

    def add_params_group_box(self):
        self.num_sigmoid_lines += 1
        params_group_box = QGroupBox()
        h_layout = QHBoxLayout()
        h_layout.addWidget(QLabel("Params #" + str(self.num_sigmoid_lines)))
        maximum = max(self.controller.scans[0].ave_smps_diameters)
        begin_rise_dp = LabelSpinbox("begin rise dp")
        begin_rise_dp.setMaximum(maximum)
        end_rise_dp = LabelSpinbox("end rise dp")
        end_rise_dp.setMaximum(maximum)
        begin_asymp_dp = LabelSpinbox("begin asymp dp")
        begin_asymp_dp.setMaximum(maximum)
        end_asymp_dp = LabelSpinbox("end asymp dp")
        end_asymp_dp.setMaximum(maximum)
        h_layout.addWidget(begin_rise_dp)
        h_layout.addWidget(end_rise_dp)
        h_layout.addWidget(begin_asymp_dp)
        h_layout.addWidget(end_asymp_dp)
        params_group_box.setLayout(h_layout)
        self.dp_widgets.append(params_group_box)
        self.params_container.addWidget(params_group_box)

    def apply(self):
        param_list = []
        for a_param_set in self.dp_widgets:
            begin_rise = a_param_set[0].content_box.value()
            end_rise = a_param_set[1].content_box.value()
            begin_asymp = a_param_set[2].content_box.value()
            end_asymp = a_param_set[3].content_box.value()
            param_list.append([begin_rise, end_rise, begin_asymp, end_asymp])
        self.controller.scans[self.controller.curr_scan_index].set_sigmoid_params(param_list)
        self.controller.switch_to_scan(self.controller.curr_scan_index)
        self.accept()


class SmoothAlgoDialog(QDialog):
    def __init__(self, controller):
        super(self.__class__, self).__init__()
        # the master layout
        form_layout = QFormLayout()
        self.scans = controller.scans
        self.controller = controller
        self.setWindowTitle("Select smoothing algorithm")
        # set up the groupbox which contains the graph of the scan
        # we assume that a scan in the middle will look very nice.
        self.curr_scan_index = len(self.scans) * 2 // 3
        while not self.scans[self.curr_scan_index].is_valid():
            self.curr_scan_index += 1
        self.a_good_scan = None
        self.original_graph = SelectSmoothingAlgorithmGraph(self)
        self.smooth_graph = SelectSmoothingAlgorithmGraph(self)
        h_layout = QHBoxLayout()
        h_layout.addWidget(self.original_graph)
        h_layout.addWidget(self.smooth_graph)
        graph_group_box = QGroupBox()
        graph_group_box.setLayout(h_layout)
        form_layout.addRow(graph_group_box)
        # set up the controls
        # self.scan_control_box = LabelSpinBox("Prev", "Next",forward=True)
        self.scan_control_box = ArrowSpinBox(forward=True)
        self.scan_control_box.setCallback(self.scan_index_changed)
        self.scan_control_box.setRange(0, len(self.scans) - 1)
        self.scan_control_box.setValue(self.curr_scan_index)
        form_layout.addRow("Select scan number", self.scan_control_box)
        # the dropdown which shows which smoothing algorithm to use
        smooth_alga_combo_box = QComboBox()
        algos = settings.smooth_algos
        for algo in algos:
            smooth_alga_combo_box.addItem(algo)
        form_layout.addRow("Select smoothing algorithm", smooth_alga_combo_box)
        # the final accept button
        button_boxes = QDialogButtonBox()
        next_button = QPushButton(self.tr("&Apply"))
        next_button.setAutoDefault(False)
        next_button.clicked.connect(self.apply)
        button_boxes.addButton(next_button, QDialogButtonBox.ApplyRole)
        form_layout.addWidget(button_boxes)
        self.setLayout(form_layout)
        # let's load up a scan
        self.process_a_scan()

    @Slot(int)
    def scan_index_changed(self, index):
        self.curr_scan_index = index
        self.process_a_scan()

    def process_a_scan(self):
        # We will deep copy it, so that we don't modify it
        # self.a_good_scan = copy.deepcopy(self.scans[self.curr_scan_index])
        self.a_good_scan = deepcopy(self.scans[self.curr_scan_index])
        # We got to produce the right amount of data first
        self.original_graph.update_graph(self.a_good_scan)
        self.smooth_graph.update_graph(self.a_good_scan, smooth_ccnc=True)

    def apply(self):
        self.accept()


class SetBaseShiftDialog(QDialog):
    def __init__(self, controller):
        super(self.__class__, self).__init__()
        # the master layout
        form_layout = QFormLayout()
        self.scans = controller.scans
        self.controller = controller
        self.curr_shift_factor = 0
        # set up the groupbox which contains the graph of the scan
        # we assume that a scan in the middle will look very nice.
        self.curr_scan_index = 0
        self.setWindowTitle("Chose a shift factor")
        while not self.scans[self.curr_scan_index].is_valid():
            self.curr_scan_index += 1
        self.a_good_scan = None
        self.shift_graph = ConcOverTimeSmoothGraph(self)
        h_layout = QHBoxLayout()
        h_layout.addWidget(self.shift_graph)
        graph_group_box = QGroupBox()
        graph_group_box.setLayout(h_layout)
        form_layout.addRow(graph_group_box)
        # set up the groupbox which contains the control for the current scan to show
        # self.scan_control_box = LabelSpinBox("Prev", "Next", forward=True)
        self.scan_control_box = ArrowSpinBox(forward=True)
        self.scan_control_box.setCallback(self.scan_index_changed)
        self.scan_control_box.setRange(0, len(self.scans) - 1)
        self.scan_control_box.setValue(self.curr_scan_index)
        form_layout.addRow("Select scan number", self.scan_control_box)
        # the area to control shifting
        self.process_a_scan()
        # self.shift_control_box = LabelSpinBox("+1S","-1S",forward=False)
        self.shift_control_box = self.scan_control_box = ArrowSpinBox(forward=True)
        self.shift_control_box.setCallback(self.shift_factor_changed)
        self.shift_control_box.setRange(-self.scans[0].duration // 2, self.scans[0].duration // 2)
        self.shift_control_box.setValue(0)
        form_layout.addRow("Select shift factor (s)", self.shift_control_box)

        # the final accept button
        button_boxes = QDialogButtonBox()
        next_button = QPushButton(self.tr("&Apply"))
        next_button.setDefault(False)
        next_button.setAutoDefault(False)
        next_button.clicked.connect(self.apply)
        button_boxes.addButton(next_button, QDialogButtonBox.ApplyRole)
        form_layout.addWidget(button_boxes)
        self.setLayout(form_layout)

    @Slot(int)
    def scan_index_changed(self, index):
        self.curr_scan_index = index
        self.process_a_scan()

    @Slot(int)
    def shift_factor_changed(self, n):
        self.curr_shift_factor = n
        self.a_good_scan.set_shift_factor(self.curr_shift_factor)
        self.a_good_scan.generate_processed_data()
        self.shift_graph.update_graph(self.a_good_scan)

    def process_a_scan(self):
        # We will deep copy it, so that we don't modify it
        self.a_good_scan = deepcopy(self.scans[self.curr_scan_index])
        self.a_good_scan.set_shift_factor(self.curr_shift_factor)
        # We got to produce the right amount of data first
        self.a_good_scan.generate_processed_data()
        self.shift_graph.update_graph(self.a_good_scan)

    def apply(self):
        self.accept()
        dialog = QMessageBox()
        dialog.setWindowTitle("Confirm Auto Alignment")
        dialog.setText("Are you sure you want to let the program to automatically align the data?")
        dialog.setIcon(QMessageBox.Question)
        dialog.addButton("Cancel", QMessageBox.RejectRole)
        confirm_button = dialog.addButton("Confirm", QMessageBox.AcceptRole)
        dialog.exec_()
        if dialog.clickedButton() == confirm_button:
            self.controller.set_base_shift_factor(self.curr_shift_factor)
            self.controller.align_smps_ccnc_data()
