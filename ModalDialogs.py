from PySide.QtGui import *
from Graphs import *
import settings
import copy
import HelperFunctions


# class InputForm(QDialog):
#     def __init__(self, mainWindow=None):
#         super(self.__class__, self).__init__()
#         self.mainWindow = mainWindow
#         self.formLayout = QFormLayout()
#         self.minDpLine = QLineEdit()
#         self.minDpAsymLine = QLineEdit()
#         self.maxDpAsymLine = QLineEdit()
#         self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
#         self.formLayout.addRow(self.tr("&MinDp"), self.minDpLine)
#         self.formLayout.addRow(self.tr("&MinDpAsym"), self.minDpAsymLine)
#         self.formLayout.addRow(self.tr("&MaxDpAsym"), self.maxDpAsymLine)
#
#         self.buttonBox.accepted.connect(self.acceptInput)
#         self.buttonBox.rejected.connect(self.reject)
#
#         self.formLayout.addWidget(self.buttonBox)
#         self.setLayout(self.formLayout)
#
#     def acceptInput(self):
#         try:
#             self.minDp = float(self.minDpLine.text())
#             self.minDpAsym = float(self.minDpAsymLine.text())
#             self.maxDpAsym = float(self.maxDpAsymLine.text())
#             self.accept()
#         except:
#             self.mainWindow.show_error_dialog("Input processed_data not valid.Please input again!")
#
#     def getData(self):
#         return self.minDp, self.minDpAsym, self.maxDpAsym
#

class SmoothAlgoDialog(QDialog):
    def __init__(self, scans, controller):
        super(self.__class__, self).__init__()
        # the master layout
        form_layout = QFormLayout()
        self.scans = scans
        self.controller = controller
        self.setWindowTitle("Select smoothing algorithm")
        # set up the groupbox which contains the graph of the scan
        # we assume that a scan in the middle will look very nice.
        self.curr_scan_index = len(scans) * 2 / 3
        while not self.scans[self.curr_scan_index].is_valid():
            self.curr_scan_index += 1
        self.a_good_scan = None
        self.original_graph = ConcOverTimeRawDataGraph(self)
        self.smooth_graph = ConcOverTimeRawDataGraph(self)
        h_layout = QHBoxLayout()
        h_layout.addWidget(self.original_graph)
        h_layout.addWidget(self.smooth_graph)
        graph_group_box = QGroupBox()
        graph_group_box.setLayout(h_layout)
        form_layout.addWidget(graph_group_box)
        # set up the groupbox which contains the control for the current scan to show
        self.curr_scan_text = QLabel("Currently showing scan " + str(self.curr_scan_index))
        self.prev_scan_button = QPushButton("Previous Scan")
        self.prev_scan_button.clicked.connect(self.prev_scan)
        self.next_scan_button = QPushButton("Next Scan")
        self.next_scan_button.clicked.connect(self.next_scan)
        h_layout = QHBoxLayout()
        h_layout.addWidget(self.curr_scan_text)
        h_layout.addWidget(self.prev_scan_button)
        h_layout.addWidget(self.next_scan_button)
        control_group_box = QGroupBox()
        control_group_box.setLayout(h_layout)
        form_layout.addWidget(control_group_box)
        # the dropdown which shows which smoothing algorithm to use
        smooth_alga_combo_box = QComboBox()
        text = QLabel("Select the smoothing algorithm to use:")
        algos = settings.smooth_algos
        for algo in algos:
            smooth_alga_combo_box.addItem(algo)
        form_layout.addWidget(text)
        form_layout.addWidget(smooth_alga_combo_box)
        # the final accept button
        button_boxes = QDialogButtonBox()
        next_button = button_boxes.addButton("Next", QDialogButtonBox.ApplyRole)
        next_button.clicked.connect(self.to_next)
        form_layout.addWidget(button_boxes)
        self.setLayout(form_layout)
        # Set the original size
        self.resize(1000, 600)
        # let's load up a scan
        self.process_a_scan()
        CONST.smooth_algos = []

    def next_scan(self):
        self.curr_scan_index = min(len(self.scans) - 1, self.curr_scan_index + 1)
        self.process_a_scan()

    def prev_scan(self):
        self.curr_scan_index = max(0, self.curr_scan_index - 1)
        self.process_a_scan()

    def process_a_scan(self):
        # We will deep copy it, so that we don't modify it
        self.a_good_scan = copy.deepcopy(self.scans[self.curr_scan_index])
        # We got to produce the right amount of data first
        self.original_graph.update_graph(self.a_good_scan)
        self.smooth_graph.update_graph(self.a_good_scan, smooth_ccnc=True)
        self.curr_scan_text.setText("Currently showing scan " + str(self.curr_scan_index))

    def to_next(self):
        self.accept()
        next_dialog = SetBaseShiftDialog(self.scans, self.controller)
        next_dialog.exec_()


class SetBaseShiftDialog(QDialog):
    def __init__(self, scans, controller):
        super(self.__class__, self).__init__()
        # the master layout
        form_layout = QFormLayout()
        self.scans = scans
        self.controller = controller
        self.curr_shift_factor = 0
        # set up the groupbox which contains the graph of the scan
        # we assume that a scan in the middle will look very nice.
        self.curr_scan_index = 0
        self.setWindowTitle("Chose base shift factor")
        while not self.scans[self.curr_scan_index].is_valid():
            self.curr_scan_index += 1
        self.a_good_scan = None
        self.original_graph = ConcOverTimeRawDataGraph(self)
        self.shift_graph = ConcOverTimeSmoothGraph(self)
        h_layout = QHBoxLayout()
        h_layout.addWidget(self.original_graph)
        h_layout.addWidget(self.shift_graph)
        graph_group_box = QGroupBox()
        graph_group_box.setLayout(h_layout)
        form_layout.addWidget(graph_group_box)
        # set up the groupbox which contains the control for the current scan to show
        self.curr_scan_text = QLabel("Currently showing scan " + str(self.curr_scan_index))
        self.prev_scan_button = QPushButton("Previous Scan")
        self.prev_scan_button.clicked.connect(self.prev_scan)
        self.next_scan_button = QPushButton("Next Scan")
        self.next_scan_button.clicked.connect(self.next_scan)
        h_layout = QHBoxLayout()
        h_layout.addWidget(self.curr_scan_text)
        h_layout.addWidget(self.prev_scan_button)
        h_layout.addWidget(self.next_scan_button)
        control_group_box = QGroupBox()
        control_group_box.setLayout(h_layout)
        form_layout.addWidget(control_group_box)
        # the area to control shifting
        self.process_a_scan()
        self.curr_shift_text = QLabel("Current shift amount is " + str(self.curr_shift_factor))
        self.shift_left_button = QPushButton("Shift Left (+1s)")
        self.shift_left_button.pressed.connect(self.shift_left)
        self.shift_left_button.setAutoRepeat(True)
        self.shift_right_button = QPushButton("shift Right (-1s)")
        self.shift_right_button.setAutoRepeat(True)
        self.shift_right_button.pressed.connect(self.shift_right)
        h_layout = QHBoxLayout()
        h_layout.addWidget(self.curr_shift_text)
        h_layout.addWidget(self.shift_left_button)
        h_layout.addWidget(self.shift_right_button)
        shift_group_box = QGroupBox()
        shift_group_box.setLayout(h_layout)
        form_layout.addWidget(shift_group_box)
        # the final accept button
        button_boxes = QDialogButtonBox()
        next_button = button_boxes.addButton("Apply", QDialogButtonBox.ApplyRole)
        next_button.clicked.connect(self.to_next)
        form_layout.addWidget(button_boxes)
        self.setLayout(form_layout)
        # Set the original size
        self.resize(1000, 600)
        # let's load up a scan

    def next_scan(self):
        self.curr_scan_index = min(len(self.scans) - 1, self.curr_scan_index + 1)
        self.process_a_scan()

    def prev_scan(self):
        self.curr_scan_index = max(0, self.curr_scan_index - 1)
        self.process_a_scan()

    def process_a_scan(self):
        # We will deep copy it, so that we don't modify it
        self.a_good_scan = copy.deepcopy(self.scans[self.curr_scan_index])
        self.a_good_scan.set_shift_factor(self.curr_shift_factor)
        # We got to produce the right amount of data first
        self.a_good_scan.generate_processed_data()
        self.original_graph.update_graph(self.a_good_scan, smooth_ccnc=True)
        self.shift_graph.update_graph(self.a_good_scan)
        self.curr_scan_text.setText("Currently showing scan " + str(self.curr_scan_index))

    def shift_left(self):
        self.curr_shift_factor += 1
        self.a_good_scan.set_shift_factor(self.curr_shift_factor)
        self.a_good_scan.generate_processed_data()
        self.shift_graph.update_graph(self.a_good_scan)
        self.curr_shift_text.setText("Current shift amount is " + str(self.a_good_scan.shift_factor))

    def shift_right(self):
        self.curr_shift_factor -= 1
        self.a_good_scan.set_shift_factor(self.curr_shift_factor)
        self.a_good_scan.generate_processed_data()
        self.shift_graph.update_graph(self.a_good_scan)
        self.curr_shift_text.setText("Current shift amount is " + str(self.a_good_scan.shift_factor))

    def to_next(self):
        self.accept()
        dialog = QMessageBox()
        dialog.setWindowTitle("Confirm Auto Alignment")
        dialog.setText("Are you sure you want to let the program to automatically align the data?")
        dialog.setIcon(QMessageBox.Question)
        cancel_button = dialog.addButton("Cancel", QMessageBox.RejectRole)
        confirm_button = dialog.addButton("Confirm", QMessageBox.AcceptRole)
        dialog.exec_()
        if dialog.clickedButton() == confirm_button:
            self.controller.set_base_shift_factor(self.curr_shift_factor)
            self.controller.align_smps_ccnc_data()
