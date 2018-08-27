from CustomWidgets import *
from datetime import date
from CustomModalDialogs import ScanInformationDialog

class DockerWidgetAlignment(QFrame):
    def __init__(self, controller):
        super(self.__class__, self).__init__()
        self.controller = controller
        # set up the layout
        form_layout = QFormLayout()
        self.setContentsMargins(-10, -10, 0, -10)
        form_layout.setContentsMargins(30, 20, 20, 0)
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Background, Qt.white)
        self.setPalette(palette)
        self.setFrameShape(QFrame.StyledPanel)
        self.setFrameShadow(QFrame.Plain)
        # The basic information of a scan
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
        # -- add total numbe rof scan
        self.num_scan = QLabel("18")
        self.num_scan.setAlignment(Qt.AlignRight)
        form_layout.addRow("Number of scan", self.num_scan)
        # -- add Counts2ConcConv
        self.counts_2_conc_conv = QLabel("0.06")
        self.counts_2_conc_conv.setAlignment(Qt.AlignRight)
        form_layout.addRow("Counts2ConcConv", self.counts_2_conc_conv)
        # -------------------
        # add a title
        form_layout.addRow(TitleHLine("Scan Information"))
        # -- add show data action
        self.show_data_button = QPushButton("Show Data")
        self.show_data_button.clicked.connect(self.show_data)
        form_layout.addRow("Data", self.show_data_button)
        # -- add the scan time
        self.scan_time = QLabel("10:26:25-10:26:25")
        self.scan_time.setAlignment(Qt.AlignRight)
        form_layout.addRow("Scan Time (h:m:s)", self.scan_time)
        # -- add the scan selector
        # self.scan_selector = LabelSpinBox("Prev","Next",forward=True)
        self.scan_selector = ArrowSpinBox(forward=True)
        self.scan_selector.setCallback(self.scan_index_changed)
        form_layout.addRow("Scan number", self.scan_selector)
        # -- add the shift selector
        self.shift_selector = ArrowSpinBox(forward=True)
        self.shift_selector.setCallback(self.shift_factor_changed)
        form_layout.addRow("Shift (s)", self.shift_selector)
        # -- add the super saturation indicator
        self.super_saturation = QLabel("0.2")
        self.super_saturation.setAlignment(Qt.AlignRight)
        form_layout.addRow("Super Saturation (%)", self.super_saturation)
        # -- add the status of the scan
        self.scan_status = QLabel("VALID")
        self.scan_status.setAlignment(Qt.AlignRight)
        form_layout.addRow("Scan status", self.scan_status)
        # add additional information
        # -- add the status of the scan
        form_layout.addRow("Additional Info", None)
        self.additional_information = QTextEdit("Welcome to Chemics!")
        self.additional_information.setReadOnly(True)
        self.additional_information.setAlignment(Qt.AlignLeft)
        form_layout.addRow(self.additional_information)
        # -- add the enable/disable button
        self.enable_disable_button = QPushButton("Disable this scan")
        self.enable_disable_button.clicked.connect(self.set_scan_enable_status)
        form_layout.addRow(self.enable_disable_button)
        # add the layout
        self.setLayout(form_layout)

    def show_data(self):
        a_scan = self.controller.scans[self.controller.curr_scan_index]
        dialog =  ScanInformationDialog(a_scan)
        dialog.exec_()

    def update_scan_info(self):
        a_scan = self.controller.scans[self.controller.curr_scan_index]
        self.counts_2_conc_conv.setText(str(self.controller.counts_to_conc_conv))
        start_time = date.strftime(a_scan.start_time, "%H:%M:%S")
        end_time = date.strftime(a_scan.end_time, "%H:%M:%S")
        scan_time = start_time + " - " + end_time
        self.scan_time.setText(scan_time)
        self.super_saturation = a_scan.processed_super_sats
        self.scan_selector.setValue(self.controller.curr_scan_index)
        self.shift_selector.setValue(a_scan.shift_factor)
        if a_scan.is_valid():
            self.scan_status.setText("VALID")
            self.scan_status.setStyleSheet("QWidget { background-color:None}")
            self.additional_information.setText("The scan shows no problem.")
            self.enable_disable_button.setText("Disable this scan")
        else:
            self.scan_status.setText("INVALID")
            self.scan_status.setStyleSheet("QWidget { color: white; background-color:red}")
            self.additional_information.setText(a_scan.decode_status_code())
            self.enable_disable_button.setText("Enable this scan")

    def update_experiment_info(self):
        self.experiment_date.setText(self.controller.experiment_date)
        self.num_scan = len(self.controller.scans) - 1
        self.shift_selector.setRange(-self.controller.scans[0].duration // 2, self.controller.scans[0].duration // 2)
        self.scan_selector.setRange(0, self.num_scan)

    def scan_index_changed(self, n):
        self.controller.set_scan_index(n)
        self.controller.switch_to_scan(n)

    def shift_factor_changed(self, n):
        curr_scan = self.controller.scans[self.controller.curr_scan_index]
        curr_scan.set_shift_factor(n)
        curr_scan.generate_processed_data()
        self.controller.switch_to_scan(self.controller.curr_scan_index)

    def set_disable_shift(self, disable = False):
        self.shift_selector.setDisabled(disable)

    def set_scan_enable_status(self):
        curr_scan = self.controller.scans[self.controller.curr_scan_index]
        if curr_scan.status == 1:
            # the user manually sets the scan to bad, so got to indicate than
            curr_scan.status = 1 - curr_scan.status
            curr_scan.status_code = 9
        else:
            # otherwise, sets it from bad to good
            curr_scan.status = 1 - curr_scan.status
        if curr_scan.is_valid():
            self.scan_status.setText("VALID")
            self.scan_status.setStyleSheet("QWidget { background-color:None}")
            self.additional_information.setText("The scan shows no problem.")
            self.enable_disable_button.setText("Disable this scan")
        else:
            self.scan_status.setText("INVALID")
            self.scan_status.setStyleSheet("QWidget { color: white; background-color:red}")
            self.additional_information.setText(curr_scan.decode_status_code())
            self.enable_disable_button.setText("Enable this scan")

class DockerSigmoidWidget(QFrame):
    def __init__(self, controller):
        super(self.__class__, self).__init__()
        # set up the layout
        form_layout = QFormLayout()
        self.setContentsMargins(-10, -10, 0, -10)
        form_layout.setContentsMargins(30, 20, 20, 0)
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Background, Qt.white)
        self.setPalette(palette)
        self.setFrameShape(QFrame.StyledPanel)
        self.setFrameShadow(QFrame.Plain)
        # next is the real data
        self.controller = controller
        self.num_sigmoid_lines = 0
        self.dp_widgets = []
        self.curr_scan_index = self.controller.curr_scan_index
        # add places to control scan and shift
        self.scan_selector = ArrowSpinBox(forward=True)
        self.scan_selector.setCallback(self.scan_index_changed)
        form_layout.addRow("Scan number", self.scan_selector)
        # -- add the status of the scan
        self.scan_status = QLabel("VALID")
        self.scan_status.setAlignment(Qt.AlignRight)
        form_layout.addRow("Scan Status", self.scan_status)
        # -- add the area to select params for sigmoid
        self.sigmoid_line_spinbox = ArrowSpinBox(forward=True)
        form_layout.addRow("Number of sigmoid lines", self.sigmoid_line_spinbox)
        self.sigmoid_line_spinbox.setCallback(self.num_params_changed)
        # add a title
        # add the items
        # add the layout
        button_boxes = QDialogButtonBox()
        apply_button = button_boxes.addButton("Apply", QDialogButtonBox.ApplyRole)
        apply_button.clicked.connect(self.apply)
        form_layout.addWidget(button_boxes)
        self.setLayout(form_layout)

    def num_params_changed(self,n):
        if n > self.num_sigmoid_lines:
            for i in range(n-self.num_sigmoid_lines):
                self.add_params_group_box()
        else:
            for i in range(self.num_sigmoid_lines - n):
                self.sub_params_group_box()

    def sub_params_group_box(self):
        self.num_sigmoid_lines = max(self.num_sigmoid_lines - 1, 0)
        self.sigmoid_line_spinbox.setValue(self.num_sigmoid_lines)
        if len(self.dp_widgets) > 0:
            to_del = self.layout().takeAt(len(self.dp_widgets)+6)
            del self.dp_widgets[-1]
            to_del.widget().deleteLater()

    def add_params_group_box(self, sigmoid_params = [0, 0, 0, 0], dp_params = [0,0,0]):
        self.num_sigmoid_lines += 1
        self.sigmoid_line_spinbox.setValue(self.num_sigmoid_lines)
        params_group_box = QGroupBox("Params #" + str(self.num_sigmoid_lines))
        maximum = max(self.controller.scans[self.curr_scan_index].ave_smps_diameters)
        begin_rise_dp = LabelDoubleSpinbox("begin rise dp")
        begin_rise_dp.setMaximum(maximum)
        begin_rise_dp.setValue(sigmoid_params[0])
        end_rise_dp = LabelDoubleSpinbox("end rise dp")
        end_rise_dp.setMaximum(maximum)
        end_rise_dp.setValue(sigmoid_params[1])
        begin_asymp_dp = LabelDoubleSpinbox("begin asymp dp")
        begin_asymp_dp.setMaximum(maximum)
        begin_asymp_dp.setValue(sigmoid_params[2])
        end_asymp_dp = LabelDoubleSpinbox("end asymp dp")
        end_asymp_dp.setValue(sigmoid_params[3])
        end_asymp_dp.setMaximum(maximum)
        v_layout = QVBoxLayout()
        v_layout.addWidget(begin_rise_dp)
        v_layout.addWidget(end_rise_dp)
        v_layout.addWidget(begin_asymp_dp)
        v_layout.addWidget(end_asymp_dp)
        h_layout = QHBoxLayout()
        dp_group_box = QGroupBox("Dp50 parameters")
        dp_50_label = QLabel("Dp50")
        dp_50_box = QLineEdit(str(dp_params[0]))
        dp_50_box.setReadOnly(True)
        another_v_layout = QVBoxLayout()
        another_v_layout.addWidget(dp_50_label)
        another_v_layout.addWidget(dp_50_box)
        h_layout.addLayout(another_v_layout)
        dp_50_wet_label = QLabel("Dp50_wet")
        dp_50_wet_box = QLineEdit(str(dp_params[1]))
        dp_50_wet_box.setReadOnly(True)
        another_v_layout = QVBoxLayout()
        another_v_layout.addWidget(dp_50_wet_label)
        another_v_layout.addWidget(dp_50_wet_box)
        h_layout.addLayout(another_v_layout)
        dp_50_20_label = QLabel("Dp50+20 wet")
        dp_50_20_box = QLineEdit(str(dp_params[2]))
        dp_50_20_box.setReadOnly(True)
        another_v_layout = QVBoxLayout()
        another_v_layout.addWidget(dp_50_20_label)
        another_v_layout.addWidget(dp_50_20_box)
        h_layout.addLayout(another_v_layout)
        dp_group_box.setLayout(h_layout)
        v_layout.addWidget(dp_group_box)
        params_group_box.setLayout(v_layout)
        self.dp_widgets.append([begin_rise_dp,end_rise_dp,begin_asymp_dp,end_asymp_dp])
        self.layout().insertRow(len(self.dp_widgets)+2, params_group_box)

    def apply(self):
        param_list = []
        for a_param_set in self.dp_widgets:
            begin_rise = a_param_set[0].content_box.value()
            end_rise = a_param_set[1].content_box.value()
            begin_asymp = a_param_set[2].content_box.value()
            end_asymp = a_param_set[3].content_box.value()
            param_list.append([begin_rise,end_rise,begin_asymp,end_asymp])
        # set the sigmoid parameters. Also immediately fit new sigmoid lines using the new params
        self.controller.scans[self.controller.curr_scan_index].set_sigmoid_params(param_list)
        self.controller.switch_to_scan(self.controller.curr_scan_index)

    def scan_index_changed(self, n):
        self.controller.set_scan_index(n)
        self.controller.switch_to_scan(n)

    def update_scan_info(self):
        num_scan = len(self.controller.scans) - 1
        self.scan_selector.setRange(0, num_scan)
        curr_scan = self.controller.scans[self.controller.curr_scan_index]
        self.scan_selector.setValue(self.controller.curr_scan_index)
        if curr_scan.is_valid():
            self.scan_status.setText("VALID")
            self.scan_status.setStyleSheet("QWidget { background-color:None}")
        else:
            self.scan_status.setText("INVALID")
            self.scan_status.setStyleSheet("QWidget { color: white; background-color:red}")
        # remove all dp
        for i in range(len(self.dp_widgets)):
            self.sub_params_group_box()
        self.num_sigmoid_lines = 0
        sigmoid_params = curr_scan.sigmoid_params
        dp_params = curr_scan.dps
        for i in range(len(sigmoid_params)):
            a_sigmoid_param = sigmoid_params[i]
            a_dp_param = [0,0,0]
            if i < len(dp_params):
                a_dp_param = dp_params[i]
            self.add_params_group_box(a_sigmoid_param, a_dp_param)


class DockerKappaWidget(QFrame):
    def __init__(self, controller, kappa_graph):
        super(self.__class__, self).__init__()
        # set up the layout
        v_layout = QVBoxLayout()
        self.setContentsMargins(-10, -10, 0, -10)
        v_layout.setContentsMargins(30, 20, 20, 0)
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Background, Qt.white)
        self.setPalette(palette)
        self.setFrameShape(QFrame.StyledPanel)
        self.setFrameShadow(QFrame.Plain)
        self.controller = controller
        self.kappa_graph = kappa_graph
        # next is the real data

        # groupbox to control showing k lines
        show_k_lines_groupbox = QGroupBox("Kappa lines")
        self.show_all_lines_radio_button = QRadioButton("Show all lines")
        self.show_all_lines_radio_button.clicked.connect(self.show_all_k_lines)
        self.show_tight_lines_radio_button =QRadioButton("Show tight lines")
        self.show_tight_lines_radio_button.clicked.connect(self.show_tight_k_lines)
        self.show_all_lines_radio_button.setChecked(True)
        h_layout = QHBoxLayout()
        h_layout.addWidget(self.show_all_lines_radio_button)
        h_layout.addWidget(self.show_tight_lines_radio_button)
        show_k_lines_groupbox.setLayout(h_layout)
        v_layout.addWidget(show_k_lines_groupbox)
        # groupbox to control showing k points
        show_k_points_groupbox = QGroupBox("Kappa values")
        self.show_all_points_radio_button = QRadioButton("Show all Ks")
        self.show_all_points_radio_button.clicked.connect(self.show_all_k_points)
        self.show_ave_points_radio_button = QRadioButton("Show averge Ks")
        self.show_ave_points_radio_button.clicked.connect(self.show_ave_k_points)
        self.show_all_points_radio_button.setChecked(True)
        h_layout = QHBoxLayout()
        h_layout.addWidget(self.show_all_points_radio_button)
        h_layout.addWidget(self.show_ave_points_radio_button)
        show_k_points_groupbox.setLayout(h_layout)
        v_layout.addWidget(show_k_points_groupbox)
        self.kappa_data_table = KappaTableWidget(self)
        v_layout.addWidget(self.kappa_data_table)
        # set layout
        self.setLayout(v_layout)

    def show_all_k_lines(self):
        self.kappa_graph.update_all_klines()

    def show_tight_k_lines(self):
        self.kappa_graph.update_tight_klines(self.controller.alpha_pinene_dict)

    def show_all_k_points(self):
        self.kappa_graph.update_all_kappa_points(self.controller.alpha_pinene_dict,
                                                 self.controller.is_valid_kappa_points)
    def show_ave_k_points(self):
        self.kappa_graph.update_average_kappa_points(self.controller.alpha_pinene_dict)

    def update_kappa_graph(self):
        if self.show_all_points_radio_button.isChecked():
            self.show_all_k_points()
        else:
            self.show_ave_k_points()

    def toggle_k_points(self, ss, dp, state):
        ss = float(ss)
        dp = float(dp)
        self.controller.set_kappa_point_state(ss,dp,state)
        self.kappa_data_table.setStatus(ss,dp,self.controller.is_valid_kappa_points[(dp,ss)])

    def update_kappa_values(self):
        for a_key in self.controller.kappa_calculate_dict.keys():
            a_scan = self.controller.kappa_calculate_dict[a_key]
            for aSS in a_scan:
                ss = a_key
                dp_50 = aSS[0]
                app_k = aSS[1]
                self.kappa_data_table.addRow(ss,dp_50,app_k,self.controller.is_valid_kappa_points[(dp_50,ss)])