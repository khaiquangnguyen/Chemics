from MainView import *
import matplotlib.pyplot as plt
import matplotlib
import scipy.optimize as opt
import scipy.stats
import scipy.signal
from Exceptions import *
import gc
from HelperFunctions import *
import FastDpCalculator
import settings as CONST
import matplotlib.patches as mpatches
import os
from PySide.QtCore import SIGNAL
from PySide.QtCore import *

from OneScan import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
matplotlib.style.use('ggplot')



class Controller:
    """
    The controller class. This class is to handle all of the functionalities of the program.
    """

    def __init__(self, view):
        # the view of the program
        self.view = view
        # list of all scans
        self.scans = []
        # counts2concconv constant
        self.counts_to_conc_conv = 0
        # ThreadPool for our little progress bar
        self.thread_pool = QThreadPool()
        self.cancelling_progress_bar = False
        # called "diameter midpoint" in the SMPS file
        self.diameter_midpoints = []
        # smps and ccnc files
        self.data_files = None
        # ccnc = Cloud Condensation Nuclei Counter
        self.ccnc_data = None
        # smps = Scanning Mobility Particle Sizer
        self.smps_data = None
        # smoothing method
        self.smooth_method = CONST.smooth_algos[0]
        # base shift factor. Very useful for auto alignment
        self.base_shift_factor = 0
        # other necessaries variables for the program. they should explain themselves easily
        # however, understanding the naming of the variables requires a decent understanding of the process
        # behind them. Therefor, I suggest getting a decent understanding of how the underlying process works
        # before diving into the code
        self.current_scan = 0
        self.number_of_scan = 0
        self.scan_start_time_list = None
        self.scan_end_time_list = None
        self.experiment_date = None
        self.scan_duration = 0
        self.ave_smps_diameters = []
        self.particle_diameter_list = None
        self.normalized_concentration_list = None
        self.super_saturation_list = []
        self.is_usable_for_sigmoid_fit_list = []
        self.min_pos_SMPS_list = []
        self.min_pos_CCNC_list = []
        self.shift_factor_list = []
        self.finish_scan_alignment_and_auto_sig_fit = False
        self.flow_rate = 0.3
        self.ccn_list = None
        self.cn_list = None
        self.ccn_fixed_list = None
        self.cn_fixed_list = None
        self.g_cn_list = None
        self.g_ccn_list = None
        self.ccn_normalized_list = []
        self.ccn_cn_sim_list = []

        # variables for sigmoid fit
        self.b = 0
        self.d = 0
        self.c = 0
        self.min_ccn = 4
        self.min_dp = 0
        self.max_dp = 0
        self.max_dp_asym = 0
        self.min_dp_asym = 0
        self.dp50 = 0
        self.super_saturation_rate = 0
        self.dp50_less_count = 0
        self.dp50_more_count = 0
        self.dp50_wet = 0
        self.dp50_plus_20 = 0
        self.min_dp_list = []
        self.min_dp_asym_list = []
        self.max_dp_asym_list = []
        self.b_list = []
        self.d_list = []
        self.c_list = []
        self.dp50_list = []
        self.dp50_wet_list = []
        self.dp50_plus_20_list = []
        self.unfinished_sigmoid_fit_scans_list = []
        self.ccnc_sig_list = []
        self.ccnc_sig_list_list = []

        # for ccn/cn graph and sigmoid fit graph
        self.ccn_cn_ratio_figure = None
        self.ccn_cn_ratio_ax = None
        self.sigmoid_fit_figure = None
        self.sigmoid_fit_ax = None
        self.ccn_cn_points = None
        self.normalized_concentration_points = None
        self.ccn_cn_ratio_corrected_points = None
        self.sigmoid_line = None

        # for all scans summary graph
        self.all_scans_alignment_figure = None
        self.all_scans_alignment_ax = None
        self.all_scans_alignment_bars = []
        self.all_scans_alignment_visited_list = []

        # variables for calculating kappa.
        self.sigma = 0.072
        self.temp = 298.15
        self.aParam = 0.00000869251 * self.sigma / self.temp
        self.dd = 280
        self.iKappa = 0.00567
        self.sc = 0
        self.asc = 0
        self.dd2 = 100
        self.iKappa2 = 0.6
        self.solubility = 0.03
        self.sc2 = 0
        self.appKappa = 0
        self.anaKappa = 0
        self.trueSC = 0
        self.kappa_excel = None
        self.kappa_points_data_list = []  # format of data is (dp,ss)
        self.klines_data = None

        # for drawing kappa graph
        self.kappa_figure = None
        self.kappa_ax = None
        self.is_show_all_k_lines = False
        self.k_lines_list = []
        self.kappa_points = None
        self.max_kappa = None
        self.min_kappa = None
        self.is_show_all_k_points = True
        self.kappa_calculate_dict = {}
        self.alpha_pinene_dict = {}
        self.scan_position_in_data = []
        self.is_usable_for_kappa_cal_list = []
        # the point of graph
        self.current_kappa_point = None
        # the index of the kappa point currently selected
        self.current_kappa_point_index = None
        # whether a kappa point is included in calculating the average k
        self.kappa_points_is_included_list = {}  # format of key is (dp,ss)

    ##############################################
    #
    # PARSE FILES AND MATCH SMPS AND CCNC DATA
    #
    ##############################################

    def parse_files(self):
        smps_txt_files = []
        ccnc_csv_files = []
        # acquire the smps and ccnc files from the input files
        for a_file in self.data_files:
            if a_file.lower().endswith('.txt'):
                smps_txt_files.append(a_file)
            elif a_file.lower().endswith('.csv'):
                ccnc_csv_files.append(a_file)
        # if (len(self.smps_txt_files) == 0) or (len(self.ccnc_csv_files) == 0):
        #     raise FileNotFoundError()
        # else:
        smps_txt_files = [str(x) for x in smps_txt_files]
        smps_txt_files = smps_txt_files[0]
        ccnc_csv_files = [str(x) for x in ccnc_csv_files]
        self.ccnc_data = process_csv_files(ccnc_csv_files)
        self.smps_data = process_text_files(smps_txt_files)

    def get_scan_timestamps(self):
        scan_down_time = 0
        scan_up_time = 0
        for i in range(len(self.smps_data)):
            if ''.join(self.smps_data[i][0].split()).lower() == "scanuptime(s)":
                scan_up_time = int(self.smps_data[i][1])
                scan_down_time = int(self.smps_data[i + 1][1])
                break
        scan_duration = scan_up_time + scan_down_time
        print ("Scan duration:", scan_duration)
        scan_start_times = self.smps_data[0]
        # generate scan end times from scan start times. Also, create a scan object for each scan
        for i in range(len(scan_start_times)):
            a_scan = OneScan()
            self.scans.append(a_scan)
            start_time = datetime.strptime(scan_start_times[i], "%H:%M:%S")
            end_time = start_time + timedelta(seconds=scan_duration)
            a_scan.set_start_time(start_time)
            a_scan.set_end_time(end_time)
            a_scan.set_up_time(scan_up_time)
            a_scan.set_down_time(scan_down_time)
            a_scan.set_duration(scan_duration)
            # todo: think about what we can do to resolve these sort of conflict

    def get_normalized_concentration(self):
        """
        Get the normalized flow_rate data. This information is taken directly from the SMPS (.txt) file
        Normalized flow_rate is also called dN/dLogDp, which is the notation we use for the graph.
        For better understanding, you can read the following document.
        http://www.tsi.com/uploadedFiles/Product_Information/Literature/Application_Notes/PR-001-RevA_Aerosol
        -Statistics-AppNote.pdf
        :return:
        """
        start_line_index = 1
        end_line_index = 0
        for i in range(1, len(self.smps_data)):
            if re.search('[aParam-zA-Z]', self.smps_data[i][0]):
                end_line_index = i
                break
        for i in range(start_line_index, end_line_index):
            self.diameter_midpoints.append(self.smps_data[i][0])
            for j in range(len(self.scans)):
                a_scan = self.scans[j]
                a_scan.add_to_raw_normalized_concs(self.smps_data[i][j + 1])

    def get_smps_counts(self):
        start_line_index = 0
        end_line_index = len(self.smps_data) - 1
        for i in range(3, len(self.smps_data)):
            if re.search('[aParam-zA-Z]', self.smps_data[i][0]):
                for j in range(i + 1, len(self.smps_data)):
                    if not re.search('[aParam-zA-Z]', self.smps_data[j][0]):
                        start_line_index = j
                        break
                break
        target_time = 1
        curr_line_index = start_line_index
        count_by_scans = [0] * len(self.scans)
        sum_diameter = 0
        while True:
            curr_time = float(self.smps_data[curr_line_index][0])
            for j in range(0, len(self.scans)):
                diameter = float(self.smps_data[curr_line_index][j * 2 + 1])
                count = int(self.smps_data[curr_line_index][j * 2 + 2])
                sum_diameter += diameter * count
                count_by_scans[j] += count
            if compare_float(curr_time, target_time) or curr_line_index == end_line_index:
                target_time += 1
                if sum(count_by_scans) == 0:
                    ave_diameter = float(self.smps_data[curr_line_index][1])
                else:
                    ave_diameter = safe_div(sum_diameter, sum(count_by_scans))
                for j in range(0, len(self.scans)):
                    self.scans[j].add_to_raw_smps_counts(count_by_scans[j])
                    self.scans[j].add_to_ave_smps_diameters(ave_diameter)
                count_by_scans = [0] * len(self.scans)
                sum_diameter = 0
            curr_line_index += 1
            if curr_line_index >= end_line_index:
                break

    def get_ccnc_counts(self):
        self.ccnc_data = self.ccnc_data[1]
        first_bin_column_index = 25
        bin_sizes = create_size_list()
        # Get the first position of CCNC count in the ccnc file
        curr_scan = 0
        curr_scan_start_time = self.scans[curr_scan].start_time
        # the index at which ccnc data is in sync with smps data
        ccnc_index = 0
        while True:
            curr_ccnc_time = datetime.strptime(self.ccnc_data[ccnc_index][0], "%H:%M:%S")
            if curr_ccnc_time > curr_scan_start_time:
                self.scans[curr_scan].set_status(0)
                self.scans[curr_scan].set_status_code(1)
                curr_scan += 1
                curr_scan_start_time = self.scans[curr_scan].start_time
            elif curr_ccnc_time < curr_scan_start_time:
                ccnc_index += 1
            else:  # the current ccnc_index is where ccnc starts being in sync with smps
                break
        finish_scanning_ccnc_data = False
        while not finish_scanning_ccnc_data:
            finish_scanning_ccnc_data = False
            a_scan = self.scans[curr_scan]
            duration = a_scan.duration
            a_scan.set_index_in_ccnc_data(ccnc_index)
            # we do one thing at a time
            for i in range(duration + duration // 4):
                curr_ccnc_index = ccnc_index + i
                # if we reach out of ccnc data bound
                if curr_ccnc_index >= len(self.ccnc_data):
                    # stop scanning ccnc data
                    finish_scanning_ccnc_data = True
                    # if we did not collect enough data for a scan, then set its status to 0
                    if i < duration:
                        a_scan.set_status(0)
                        a_scan.set_status_code(1)
                    break
                # collect a bunch of data from ccnc file
                super_saturation = self.ccnc_data[curr_ccnc_index][1]
                a_scan.add_to_raw_super_sats(super_saturation)
                ccnc_count = self.ccnc_data[curr_ccnc_index][-3]
                a_scan.add_to_raw_ccnc_counts(ccnc_count)
                total_count = 0
                total_size = 0
                for j in range(len(bin_sizes)):
                    count_in_bin = int(float(self.ccnc_data[curr_ccnc_index][first_bin_column_index + j]))
                    total_size += bin_sizes[j] * count_in_bin
                    total_count += count_in_bin
                ave_size = safe_div(total_size, total_count)
                a_scan.add_to_raw_ave_ccnc_sizes(ave_size)
                a_scan.add_to_raw_t1s_t2s_t3s(self.ccnc_data[curr_ccnc_index][5], self.ccnc_data[curr_ccnc_index][7],
                                              self.ccnc_data[curr_ccnc_index][9])
            curr_scan += 1
            # if we run of out scans to compare with ccnc data, stop scanning ccnc data
            if curr_scan >= len(self.scans):
                break
            # find the next ccnc_index
            # we got to based on the start time, since the duration values are always off
            next_scan_start_time = self.scans[curr_scan].start_time
            while True:
                curr_ccnc_time = datetime.strptime(self.ccnc_data[ccnc_index][0], "%H:%M:%S")
                if curr_ccnc_time < next_scan_start_time:
                    ccnc_index += 1
                    # if we reach out of ccnc data bound
                    if ccnc_index >= len(self.ccnc_data):
                        # stop scanning ccnc data
                        finish_scanning_ccnc_data = True
                        break
                else:
                    break

        print()

    def do_basic_trans(self):
        for i in range(len(self.scans)):
            self.scans[i].do_basic_trans()

    def pre_align_sanity_check(self):
        # Perform self test for each scan
        for i in range(len(self.scans)):
            self.scans[i].pre_align_self_test()
        # Cross validation. Basically compare the distribution of a scan with the next two
        # We know that only the first few distributions have weird data, so once it becomes right, we stop
        for i in range(len(self.scans) - 2):
            # if the scan is invalid, we skip to the next one
            if not self.scans[i].is_valid():
                continue
            # if a scan has a different dist from the next one or the one after that
            if not self.scans[i].compare_smps(self.scans[i + 1]) or not self.scans[i].compare_smps(self.scans[i + 2]):
                self.scans[i].set_status(0)
                self.scans[i].set_status_code(3)
            else:
                break

    def post_align_sanity_check(self):
        # Perform self-test for each scan
        for i in range(len(self.scans)):
            self.scans[i].post_align_self_test()
        # Check if we can make any reference about that data.
        # Detect any outliers in shift factors
        shift_factors = []
        for i in range(len(self.scans)):
            shift_factors.append(self.scans[i].shift_factor)
        print shift_factors
        # todo: outliers detection
        # todo: fix all wrong shift factors
        # basically, the idea is that a shift factors can't be too different from its neighbors
        # therefore, a radical difference indiciates that the shift factor is wrong
        # we simply fix it by its neightbors. It will likely be correct, unless the data is wrong.

    def align_smps_ccnc_data(self):
        def align_smps_thread_sequence(progress_start,progress_update):
            shift_factor = self.base_shift_factor
            progress_start.emit("Aligning SMPS and CCNC data...")
            for i in range(len(self.scans)):
                progress_update.emit(100*(i+1)//len(self.scans))
                shift_factor = self.scans[i].align_smps_ccnc_data(shift_factor)
                self.scans[i].generate_processed_data()
                self.view.update_alignment_graphs(self.scans[i])
            self.view.update_alignment_graphs(self.scans[0])
            self.post_align_sanity_check()
        started_fn = self.view.init_progress_bar
        finish_fn = self.view.close_progress_bar
        progress_fn = self.view.update_progress_bar
        fns_to_run = align_smps_thread_sequence
        worker = run_new_thread(fns_to_run, started_fn, finish_fn, progress_fn)
        self.thread_pool.start(worker)


    ##############################################
    #
    # CORRECT CHARGES AND FIT SIGMOID
    #
    ##############################################

    def get_parameters_for_sigmoid_fit(self, min_dp=None, min_dp_asym=None, max_dp_asym=None):
        """
        Automatically generate parameters to fit sigmoid lines into the data. If the paramaters are declared
        then the process will work with the parameter instead of automatically deciding the parameters.
        :param min_dp:
        :param min_dp_asym:
        :param max_dp_asym:
        :return:
        """
        try:
            if min_dp and min_dp_asym and max_dp_asym:
                self.min_dp = min_dp
                self.min_dp_asym = min_dp_asym
                self.max_dp_asym = max_dp_asym
            else:
                asymList = get_asym_list(self.particle_diameter_list, self.ccnc_sig_list)
                mpIndex = 0

                # Find the index of midPoint
                checkLength = self.scan_duration / 20
                for i in range(checkLength, len(self.ccnc_sig_list) - checkLength):
                    isMid = False
                    if self.ccnc_sig_list[i] > 0.5:
                        isMid = True
                        # check if the previous 5 numbers are smaller and the next 5 numbers are bigger
                        for j in range(1, checkLength):
                            if self.ccnc_sig_list[i + j] < self.ccnc_sig_list[i]:
                                isMid = False
                                break
                    if isMid:
                        mpIndex = i
                        break

                minDpAsymPos = 0
                currMax = 0
                # Get minDP
                for i in range(mpIndex, 1, -1):
                    if self.ccnc_sig_list[i] < 0.1:
                        self.min_dp = self.particle_diameter_list[i]
                        break

                # Get minDpAsym
                for i in range(mpIndex, mpIndex + self.scan_duration / 10):
                    if i < len(self.particle_diameter_list):
                        if self.ccnc_sig_list[i] > currMax:
                            minDpAsymPos = i
                            currMax = self.ccnc_sig_list[i]
                            self.min_dp_asym = self.particle_diameter_list[i]

                maxCCN = max(self.ccnc_sig_list)
                if maxCCN <= 1.2:
                    stableThreshold = 0.075
                else:
                    stableThreshold = 0.25
                # determine maxDplen(self.ccncSigList)
                for i in range(minDpAsymPos + 5, len(self.ccnc_sig_list)):
                    if self.ccnc_sig_list[i] > 1.3:
                        self.max_dp_asym = self.particle_diameter_list[i]
                        break
                    elif abs(self.ccnc_sig_list[i] - self.ccnc_sig_list[i - 1]) > 2 * stableThreshold:
                        self.max_dp_asym = self.particle_diameter_list[i - 1]
                        break

            self.max_dp = self.max_dp_asym
            # Get the processed_data
            asymsList = []
            for i in range(len(self.particle_diameter_list)):
                if self.min_dp_asym < self.particle_diameter_list[i] < self.max_dp_asym:
                    asymsList.append(self.ccnc_sig_list[i])
                else:
                    asymsList.append(0)
            # Calcualte constants
            self.b = get_ave_none_zero(asymsList)
        except:
            if self.current_scan not in self.unfinished_sigmoid_fit_scans_list:
                self.unfinished_sigmoid_fit_scans_list.append(self.current_scan)
            raise SigmoidFitGetParameterError()

    def fit_sigmoid_line(self):
        """
        As the name suggest, using the parameters to fit a sigmoid line into the data
        I use the curve_fit function provided by scipy.optimize library to run this
        :return:
        """
        try:
            xList = []
            yList = []
            for i in range(len(self.particle_diameter_list)):
                if self.min_dp < self.particle_diameter_list[i] < self.max_dp:
                    xList.append(self.particle_diameter_list[i])
                    yList.append(self.ccnc_sig_list[i])
            initGuess = [30, -10]
            for i in range(len(yList)):
                if math.isnan(yList[i]) or math.isinf(yList[i]):
                    yList[i] = 0
            xList = numpy.asarray(xList)
            yList = numpy.asarray(yList)
            result = opt.curve_fit(f, xList, yList, bounds=([self.min_dp, -200], [self.max_dp, -1]), method="trf")
            self.d = result[0][0]
            self.c = result[0][1]
            self.ccn_cn_sim_list = [0]
            for i in range(1, len(self.particle_diameter_list)):
                if self.particle_diameter_list[i] > self.d:
                    self.dp50_wet = self.drop_size_list[i - 1]
                    break
            for i in range(1, len(self.particle_diameter_list)):
                if self.particle_diameter_list[i] > (self.d + 20):
                    self.dp50_plus_20 = self.drop_size_list[i - 1]
                    break
            self.is_usable_for_kappa_cal_list[self.current_scan] = True
        except:
            if self.current_scan not in self.unfinished_sigmoid_fit_scans_list:
                self.unfinished_sigmoid_fit_scans_list.append(self.current_scan)
            raise SigmoidFitLineFitError()

    def refitting_sigmoid_line(self, min_dry_diameter, min_dry_diameter_asymptote, max_dry_diameter_asymptote):
        """
        As the name suggest, this function is for the users to refit a failed sigmoid line.
        :param min_dry_diameter:
        :param min_dry_diameter_asymptote:
        :param max_dry_diameter_asymptote:
        :return:
        """
        self.move_progress_bar_forward("Refitting sigmoid line to scan #" + str(self.current_scan + 1),
                                       max_value=2)
        try:
            if not self.is_usable_for_sigmoid_fit_list[self.current_scan]:
                self.view.show_error_dialog("Can't fit sigmoid line to current scan!")
                raise SigmoidFitLineFitError()
            self.prepare_scan_data()
            self.move_progress_bar_forward()
            self.min_dp = min_dry_diameter
            self.min_dp_asym = min_dry_diameter_asymptote
            self.max_dp_asym = max_dry_diameter_asymptote
            self.max_dp = self.max_dp_asym
            asymsList = []
            for i in range(len(self.particle_diameter_list)):
                if self.min_dp_asym < self.particle_diameter_list[i] < self.max_dp_asym:
                    asymsList.append(self.ccnc_sig_list[i])
                else:
                    asymsList.append(0)
            self.b = get_ave_none_zero(asymsList)
            self.fit_sigmoid_line()
            self.move_progress_bar_forward()
        except:
            self.is_usable_for_kappa_cal_list[self.current_scan] = False
            self.b_list[self.current_scan] = 0
            self.d_list[self.current_scan] = 0
            self.c_list[self.current_scan] = 0
            self.dp50_list[self.current_scan] = (0, 0)
            self.dp50_wet_list[self.current_scan] = 0
            self.dp50_plus_20_list[self.current_scan] = 0
            self.min_dp_list[self.current_scan] = 0
            self.min_dp_asym_list[self.current_scan] = 0
            self.max_dp_asym_list[self.current_scan] = 0
        else:
            self.is_usable_for_kappa_cal_list[self.current_scan] = True
            self.min_dp_list[self.current_scan] = self.min_dp
            self.min_dp_asym_list[self.current_scan] = self.min_dp_asym
            self.max_dp_asym_list[self.current_scan] = self.max_dp_asym
            self.b_list[self.current_scan] = self.b
            self.d_list[self.current_scan] = self.d
            self.c_list[self.current_scan] = self.c
            self.dp50_list[self.current_scan] = (self.d, self.super_saturation_list[self.current_scan])
            self.dp50_wet_list[self.current_scan] = self.dp50_wet
            self.dp50_plus_20_list[self.current_scan] = self.dp50_plus_20
        finally:
            self.move_progress_bar_forward(complete=True)
            # remove current scan from the list of unfinished scan
            for i in range(len(self.unfinished_sigmoid_fit_scans_list)):
                if self.unfinished_sigmoid_fit_scans_list[i] == self.current_scan:
                    self.unfinished_sigmoid_fit_scans_list = self.unfinished_sigmoid_fit_scans_list[:i] \
                                                             + self.unfinished_sigmoid_fit_scans_list[i + 1:]
                    break
            self.update_view()

    def correct_charges_and_fit_sigmoid_one_scan(self):
        """
        the control center to correct charges and fit sigmoid line to a single scan.

        """
        try:
            if not self.is_usable_for_sigmoid_fit_list[self.current_scan]:
                raise SigmoidFitCorrectChargesError  # this scan is not usable, so we raise a serious error
            self.ccnc_sig_list = []
            self.prepare_scan_data()
            self.correct_charges()
            self.get_parameters_for_sigmoid_fit()
            self.fit_sigmoid_line()
        except Exception as e:
            self.is_usable_for_kappa_cal_list[self.current_scan] = False
            self.b_list.append(0)
            self.d_list.append(0)
            self.c_list.append(0)
            self.dp50_list.append((0, 0))
            self.dp50_wet_list.append(0)
            self.dp50_plus_20_list.append(0)
            self.min_dp_list.append(0)
            self.min_dp_asym_list.append(0)
            self.max_dp_asym_list.append(0)
            self.ccnc_sig_list_list.append(self.ccnc_sig_list)
        else:
            self.is_usable_for_kappa_cal_list[self.current_scan] = True
            self.min_dp_list.append(self.min_dp)
            self.min_dp_asym_list.append(self.min_dp_asym)
            self.max_dp_asym_list.append(self.max_dp_asym)
            self.b_list.append(self.b)
            self.d_list.append(self.d)
            self.c_list.append(self.c)
            self.dp50_list.append((self.d, self.super_saturation_list[self.current_scan]))
            self.dp50_wet_list.append(self.dp50_wet)
            self.dp50_plus_20_list.append(self.dp50_plus_20)
            self.ccnc_sig_list_list.append(self.ccnc_sig_list)

    def correct_charges_and_fit_sigmoid_all_scans(self):
        """
        Perform charge correction and sigmoid fit for all scans.
        :return:
        """
        try:
            self.move_progress_bar_forward(max_value=self.number_of_scan)
            for i in range(0, self.number_of_scan):
                self.current_scan = i
                self.move_progress_bar_forward("Correcting charges and fitting sigmoid for scan # " + str(i + 1),
                                               value=0)
                self.correct_charges_and_fit_sigmoid_one_scan()
            self.current_scan = -1
            self.move_progress_bar_forward(complete=True)
            self.finish_scan_alignment_and_auto_sig_fit = True
            self.switch_to_scan(0)
        except ProgressBarInterruptException:
            self.view.show_error_dialog("The sigmoid fitting process is cancelled!")
            self.finish_scan_alignment_and_auto_sig_fit = False
            self.parse_and_match_smps_ccnc_data()

    # -------------------------graphs-----------------------------

    def draw_complete_sigmoid_graph(self, new_figure=None):
        """
        Make complete graph of the dry diameter after optimization and sigmodal fit
        """
        # if there is no ax, make a new one
        if self.sigmoid_fit_ax is None:
            figure, ax = plt.subplots(facecolor=settings.GRAPH_BACKGROUND_COLOR)
            ax.axes.set_frame_on(False)
            ax.grid(color=GRID_COLOR)
            ax.axhline(0, color=AX_LINE_COLOR, linewidth=2)
            ax.axvline(0, color=AX_LINE_COLOR, linewidth=4)
            ax.axhline(1, color=str(float(AX_LINE_COLOR) + 0.1), linewidth=2, linestyle='dashed')
            ax.tick_params(color=AX_TICK_COLOR, which='both', labelcolor=AX_TICK_COLOR, labelsize=AX_TICK_SIZE)
            ax.set_xlabel("Dry diameter(nm)", color=LABEL_COLOR, size=LABEL_SIZE)
            ax.set_ylabel("CCNC/SMPS ratio and dN/dLogDp", color=LABEL_COLOR, size=LABEL_SIZE)
            ax.set_title("CCNC/SMPS over Dry Diameter and Sigmoid Line", color=TITLE_COLOR, size=TITLE_SIZE)
            yLim = min(2, max(self.ccn_cn_ratio_list)) + 0.2
            ax.axes.set_ylim([-0.1, yLim])
            self.ccn_cn_ratio_points, = ax.plot(self.particle_diameter_list, self.ccn_cn_ratio_list, 'o',
                                                color=CCNC_SMPS_POINT_COLOR,
                                                mew=0.5, ms=9, label="CCNC/SMPS")
            if self.is_usable_for_sigmoid_fit_list[self.current_scan]:
                self.normalized_concentration_points, = plt.plot(self.diameter_midpoint_list, self.ccn_normalized_list,
                                                                 linewidth=4,
                                                                 color=NORMALIZED_CONCENTRATION_POINT_COLOR,
                                                                 label="dN/dLogDp")

                self.ccn_cn_ratio_corrected_points, = ax.plot(self.particle_diameter_list, self.ccnc_sig_list, 'o',
                                                              color=CCNC_SMPS_RATIO_CORRECTED_POINT_COLOR, mew=0.5,
                                                              mec="#0D47A1",
                                                              ms=9, label="CCNC/SMPS (Corrected)")
            if self.is_usable_for_kappa_cal_list[self.current_scan]:
                self.sigmoid_line, = ax.plot(self.particle_diameter_list, self.ccn_cn_sim_list, linewidth=5,
                                             color=SIGMOID_LINE_COLOR, label="Sigmodal Fit")
            handles, labels = ax.get_legend_handles_labels()
            legend = ax.legend(handles, labels, loc=5, fontsize='small')
            legend.get_frame().set_alpha(0.3)
            self.sigmoid_fit_figure = figure
            self.sigmoid_fit_ax = ax
        # otherwise, use the old one
        else:
            yLim = min(2, max(self.ccn_cn_ratio_list)) + 0.2
            self.sigmoid_fit_ax.axes.set_ylim([-0.1, yLim])
            self.ccn_cn_ratio_points.set_xdata(self.particle_diameter_list)
            self.ccn_cn_ratio_points.set_ydata(self.ccn_cn_ratio_list)
            if self.is_usable_for_sigmoid_fit_list[self.current_scan]:
                if self.normalized_concentration_points is None:
                    self.normalized_concentration_points, = plt.plot(self.diameter_midpoint_list,
                                                                     self.ccn_normalized_list,
                                                                     linewidth=4,
                                                                     color=NORMALIZED_CONCENTRATION_POINT_COLOR,
                                                                     label="dN/dLogDp")
                else:
                    self.normalized_concentration_points.set_xdata(self.diameter_midpoint_list)
                    self.normalized_concentration_points.set_ydata(self.ccn_normalized_list)
                if self.ccn_cn_ratio_corrected_points is None:
                    self.ccn_cn_ratio_corrected_points, = self.sigmoid_fit_ax.plot(self.particle_diameter_list,
                                                                                   self.ccnc_sig_list, 'o',
                                                                                   color=CCNC_SMPS_RATIO_CORRECTED_POINT_COLOR,
                                                                                   mew=0.5,
                                                                                   mec="#0D47A1",
                                                                                   ms=9, label="CCNC/SMPS (Corrected)")
                else:
                    self.ccn_cn_ratio_corrected_points.set_xdata(self.particle_diameter_list)
                    self.ccn_cn_ratio_corrected_points.set_ydata(self.ccnc_sig_list)
            else:
                if self.normalized_concentration_points is not None:
                    self.normalized_concentration_points.set_xdata([])
                    self.normalized_concentration_points.set_ydata([])
                if self.ccn_cn_ratio_corrected_points is not None:
                    self.ccn_cn_ratio_corrected_points.set_xdata([])
                    self.ccn_cn_ratio_corrected_points.set_ydata([])
            if self.is_usable_for_kappa_cal_list[self.current_scan]:
                if self.sigmoid_line is None:
                    self.sigmoid_line, = self.sigmoid_fit_ax.plot(self.particle_diameter_list, self.ccn_cn_sim_list,
                                                                  linewidth=5, color='#EF5350', label="Sigmodal Fit")
                else:
                    self.sigmoid_line.set_xdata(self.particle_diameter_list)
                    self.sigmoid_line.set_ydata(self.ccn_cn_sim_list)
            else:
                if self.sigmoid_line is not None:
                    self.sigmoid_line.set_xdata([])
                    self.sigmoid_line.set_ydata([])
            handles, labels = self.sigmoid_fit_ax.get_legend_handles_labels()
            legend = self.sigmoid_fit_ax.legend(handles, labels, loc=5, fontsize='small')
            legend.get_frame().set_alpha(0.3)

        self.view.update_alignment_and_sigmoid_fit_figures(None, self.sigmoid_fit_figure)

    def draw_all_scans_alignment_summary_graph(self):
        """
        Make a summary graph of all the runs. Used color to show the status of each scan
        Red represents fail scans.
        Blue represents good scans.
        White represents undecided scans.
        Green (or yellow? I think I am color blind) represents visited scan
        :return:
        """
        if self.all_scans_alignment_ax is None:
            for i in range(self.number_of_scan):
                self.all_scans_alignment_visited_list.append(False)

            figure, ax = plt.subplots(facecolor=settings.GRAPH_BACKGROUND_COLOR)
            ax.axes.set_frame_on(False)
            ax.grid(False)
            ax.axhline(0, color=AX_LINE_COLOR, linewidth=4)
            ax.axvline(0, color=AX_LINE_COLOR, linewidth=4)
            index = numpy.arange(1, self.number_of_scan + 1)
            ax.set_xticks(index)
            ax.set_xticklabels(index, color=AX_TICK_COLOR, size=AX_TICK_SIZE)
            ax.set_yticks(range(0, max(self.shift_factor_list)))
            ax.set_yticklabels(range(0, max(self.shift_factor_list)), color=AX_TICK_COLOR, size=AX_TICK_SIZE)
            ax.set_xlabel("Scan #", color=LABEL_COLOR, size=LABEL_SIZE)
            ax.set_ylabel("CCNC Shift (s)", color=LABEL_COLOR, size=LABEL_SIZE)
            ax.set_title("All Scans Shift", color=TITLE_COLOR, size=TITLE_SIZE)
            figure.canvas.mpl_connect('pick_event', self.on_pick)
            self.all_scans_alignment_bars = ax.bar(range(1, self.number_of_scan + 1), self.shift_factor_list,
                                                   color=SCAN_SUMMARY_USABLE_SCAN_COLOR, picker=True, align='center')
            for i in range(len(self.is_usable_for_sigmoid_fit_list)):
                if not self.is_usable_for_sigmoid_fit_list[i] or not self.is_usable_for_kappa_cal_list[i]:
                    self.all_scans_alignment_bars[i].set_facecolor(SCAN_SUMMARY_UNUSABLE_SCAN_COLOR)
            for i in self.unfinished_sigmoid_fit_scans_list:
                self.all_scans_alignment_bars[i].set_facecolor(UNDECIDED_USABILITY_COLOR)
            self.all_scans_alignment_visited_list[self.current_scan] = True
            for i in range(len(self.all_scans_alignment_visited_list)):
                if self.all_scans_alignment_visited_list[i]:
                    self.all_scans_alignment_bars[i].set_alpha(1)
                else:
                    self.all_scans_alignment_bars[i].set_alpha(0.3)
            self.all_scans_alignment_bars[self.current_scan].set_facecolor(SCAN_SUMMARY_HIGHLIGHT_COLOR)
            self.all_scans_alignment_figure = figure
            valid_patch = mpatches.Patch(color=SCAN_SUMMARY_USABLE_SCAN_COLOR, label='usable scans')
            invalid_patch = mpatches.Patch(color=SCAN_SUMMARY_UNUSABLE_SCAN_COLOR, label='unusable scans')
            undecided_patch = mpatches.Patch(color=UNDECIDED_USABILITY_COLOR, label='undecided scans')
            visted_patch = mpatches.Patch(color=SCAN_SUMMARY_HIGHLIGHT_COLOR, label='visited scans')
            legend = ax.legend(handles=[valid_patch, invalid_patch, undecided_patch, visted_patch],
                               fontsize=LEGEND_FONT_SIZE, loc=2)
            self.all_scans_alignment_ax = ax
            legend.get_frame().set_alpha(0.3)
        else:
            for i in range(len(self.is_usable_for_sigmoid_fit_list)):
                if not self.is_usable_for_sigmoid_fit_list[i] or not self.is_usable_for_kappa_cal_list[i]:
                    self.all_scans_alignment_bars[i].set_facecolor(SCAN_SUMMARY_UNUSABLE_SCAN_COLOR)
                else:
                    self.all_scans_alignment_bars[i].set_facecolor(SCAN_SUMMARY_USABLE_SCAN_COLOR)
            for i in self.unfinished_sigmoid_fit_scans_list:
                self.all_scans_alignment_bars[i].set_facecolor(UNDECIDED_USABILITY_COLOR)
            self.all_scans_alignment_visited_list[self.current_scan] = True
            for i in range(len(self.all_scans_alignment_visited_list)):
                if self.all_scans_alignment_visited_list[i]:
                    self.all_scans_alignment_bars[i].set_alpha(1)
                else:
                    self.all_scans_alignment_bars[i].set_alpha(0.3)
            self.all_scans_alignment_bars[self.current_scan].set_facecolor(SCAN_SUMMARY_HIGHLIGHT_COLOR)

        self.view.update_temp_and_min_figure(self.all_scans_alignment_figure)

    def on_pick(self, event):
        """
        When a bar on the summary graph is clicked
        :param event:
        :return:
        """
        for i in range(len(self.all_scans_alignment_bars)):
            if event.artist == self.all_scans_alignment_bars[i]:
                if i != self.current_scan:
                    self.switch_to_scan(i)

    def on_key_release(self, event):
        """
        When a key is pressed while in the summary graph
        :param event:
        :return:
        """
        key = event.key()
        # left arrow key
        if key == 16777234 or key == 16777237:
            self.switch_to_scan(max(0, self.current_scan - 1))
        # right arrow key
        elif key == 16777236 or key == 16777235:
            self.switch_to_scan(min(self.current_scan + 1, self.number_of_scan - 1))

    ##############################################
    #
    # KAPPA CALCULATION
    #
    ##############################################

    def set_parameters_for_kappa_calculation(self, sigma=None, temp=None, dd1=None, i_kappa=None, dd2=None,
                                             i_kappa_2=None, solu=None):
        """
        Update the necessary constants to calculate Kappa
        """
        if sigma:
            self.sigma = sigma
        if temp:
            self.temp = temp
        if dd1:
            self.dd = dd1
        if i_kappa:
            self.iKappa = i_kappa
        if dd2:
            self.dd2 = dd2
        if i_kappa_2:
            self.iKappa2 = i_kappa_2
        if solu:
            self.solubility = solu

    def calculate_all_kappa_values(self):
        """
        Calculate the kappa values for every scan.
        """
        if self.kappa_excel is None:
            self.kappa_excel = pandas.read_csv("kCal.csv", header=None)
        lookup = self.kappa_excel
        self.aParam = 0.00000869251 * self.sigma / self.temp
        self.asc = (exp(sqrt(4 * self.aParam ** 3 / (27 * self.iKappa * (self.dd * 0.000000001) ** 3))) - 1) * 100
        # Calculate each kappa
        firstAKappa = 0
        scCalcs = False
        for i in range(len(self.dp50_list)):
            if not self.is_usable_for_kappa_cal_list[i]:
                continue
            ss = float(self.dp50_list[i][1])
            dp50 = float(self.dp50_list[i][0])
            rowIndex = int(math.floor(dp50 - 9))
            matchRow = list(lookup.iloc[rowIndex][1:])
            valueRow = list(lookup.iloc[0][1:])
            a = get_correct_num(matchRow, ss)
            cIndex = a[1]
            a = a[0]
            if cIndex != (len(matchRow) - 1):
                b = matchRow[cIndex + 1]
                c = valueRow[cIndex]
                d = valueRow[cIndex + 1]
            else:
                c = valueRow[cIndex]
                b = 0
                d = 0
            self.appKappa = (ss - (a - (a - b) / (c - d) * c)) / ((a - b) / (c - d))
            if not scCalcs:
                # Calculate the sc Calculation
                firstAKappa = self.appKappa
                # Calcualte the first row of scCalcs
                for i in range(len(self.dp50_list)):
                    if self.dp50_list[i][0] != 0:
                        dList1 = [float(self.dp50_list[i][0]) * 0.000000001]
                        break
                for i in range(1000):
                    dList1.append(dList1[-1] * 1.005)

                # Calcualte the second row of scCalcs
                sList1 = []
                firstNum = dList1[0]
                for i in range(len(dList1)):
                    aNum = dList1[i]
                    sList1.append(
                        (aNum ** 3 - firstNum ** 3) / (aNum ** 3 - firstNum ** 3 * (1 - firstAKappa)) * math.exp(
                            self.aParam / aNum))

                # Calculate the third colum
                dList2 = [self.dd * 0.000000001]
                for i in range(1000):
                    dList2.append(dList2[-1] * 1.005)

                # Calculate the fourth column
                sList2 = []
                firstNum = dList2[0]
                for i in range(len(dList1)):
                    aNum = dList2[i]
                    sList2.append(
                        (aNum ** 3 - firstNum ** 3) / (aNum ** 3 - firstNum ** 3 * (1 - self.iKappa)) * math.exp(
                            self.aParam / aNum))

                dList3 = [self.dd2 * 0.000000001]
                for i in range(1000):
                    dList3.append(dList3[-1] * 1.005)

                kappaList = []
                firstNum = dList3[0]
                for i in range(len(dList3)):
                    aNum = dList3[i]
                    kappaList.append(min((aNum ** 3 / firstNum ** 3 - 1) * self.solubility, 1) * self.iKappa2)

                sList3 = []
                firstNum = dList3[0]
                for i in range(len(dList1)):
                    aNum = dList3[i]
                    if kappaList[i] == 0:
                        sList3.append(0)
                    else:
                        sList3.append(
                            (aNum ** 3 - firstNum ** 3) / (aNum ** 3 - firstNum ** 3 * (1 - kappaList[i])) * math.exp(
                                self.aParam / aNum))
                self.sc = (max(sList2) - 1) * 100
                self.sc2 = (max(sList3) - 1) * 100
                scCalcs = True
            self.trueSC = (max(sList1[i:]) - 1) * 100
            self.anaKappa = (4 * self.aParam ** 3) / (27 * (dp50 * 0.000000001) ** 3 * log(ss / 100 + 1) ** 2)
            kDevi = (self.appKappa - self.anaKappa) / self.appKappa * 100
            if ss in self.kappa_calculate_dict.keys():
                self.kappa_calculate_dict[ss].append([dp50, self.appKappa, self.anaKappa, kDevi, self.trueSC])
            else:
                self.kappa_calculate_dict[ss] = ([[dp50, self.appKappa, self.anaKappa, kDevi, self.trueSC]])
            self.kappa_points_is_included_list[(dp50, ss)] = True
        self.update_kappa_info_and_graph()

    def calculate_average_kappa_values(self):
        """
        Calculate the kappa values for each super saturation percentage. The values are average of all scans with the
        same super saturation
        :return:
        """
        self.alpha_pinene_dict = {}
        for a_key in self.kappa_calculate_dict.keys():
            aSSList = self.kappa_calculate_dict[a_key]
            temp_dp50List = []
            dp50List = []
            appKappaList = []
            anaKappaList = []
            meanDevList = []
            for aSS in aSSList:
                dp50List.append(aSS[0])
                if self.kappa_points_is_included_list[(aSS[0], a_key)]:
                    temp_dp50List.append(aSS[0])
                    appKappaList.append(aSS[1])
                    anaKappaList.append(aSS[2])
                    meanDevList.append(aSS[3])
            meanDp = average(temp_dp50List)
            stdDp = numpy.std(temp_dp50List)
            meanApp = average(appKappaList)
            stdApp = numpy.std(appKappaList)
            meanAna = average(anaKappaList)
            stdAna = numpy.std(anaKappaList)
            meanDev = average(meanDevList)
            devMean = (meanApp - meanAna) / meanApp * 100
            self.alpha_pinene_dict[a_key] = (
                meanDp, stdDp, meanApp, stdApp, meanAna, stdAna, meanDev, devMean, dp50List)

    def update_kappa_info_and_graph(self):
        self.calculate_average_kappa_values()
        self.draw_kappa_graph()
        self.view.update_kappa_info_and_graph()

    # --------------------- Graphs -----------------------------

    def draw_kappa_graph(self):
        """
        Probably the most complicated graph of the whole program.
        It needs to preload the kappa lines before actually drawing the kappa points.
        :return:
        """
        if self.klines_data is None:
            self.klines_data = pandas.read_csv("klines.csv", header=1)
        header = self.klines_data.columns
        diameter_list = self.klines_data[header[1]]
        kappa_list = []
        std_kappa_list = []
        self.kappa_points_data_list = []
        # update the kappa points list.
        for key in self.alpha_pinene_dict.keys():
            kappa_list.append(self.alpha_pinene_dict[key][2])
            std_kappa_list.append(self.alpha_pinene_dict[key][3])
            if self.is_show_all_k_points:
                for i in self.alpha_pinene_dict[key][-1]:
                    self.kappa_points_data_list.append((i, key))
            else:
                if not math.isnan(self.alpha_pinene_dict[key][0]):
                    self.kappa_points_data_list.append((self.alpha_pinene_dict[key][0], key))
        self.kappa_points_data_list.sort(key=lambda tup: tup[0], reverse=True)
        self.kappa_points_data_list.sort(key=lambda tup: tup[1], reverse=True)
        if self.current_kappa_point_index is None:
            self.current_kappa_point_index = len(self.kappa_points_data_list) - 1
        k_points_x_list = []
        k_points_y_list = []
        all_k_points_color_list = []
        for i in range(len(self.kappa_points_data_list)):
            # x and y pos of points
            k_points_x_list.append(self.kappa_points_data_list[i][0])
            k_points_y_list.append(self.kappa_points_data_list[i][1])
            # get the colors for the included and excluded points
            if self.is_show_all_k_points:
                ss = self.kappa_points_data_list[i][1]
                dp = self.kappa_points_data_list[i][0]
                if self.kappa_points_is_included_list[(dp, ss)]:
                    all_k_points_color_list.append(KAPPA_USABLE_POINT_COLOR)
                else:
                    all_k_points_color_list.append(KAPPA_UNUSABLE_POINT_COLOR)
            else:
                all_k_points_color_list.append(KAPPA_USABLE_POINT_COLOR)
        if self.is_show_all_k_lines:
            kappa_start_pos = 2
            kappa_end_pos = len(header)
        else:
            # Get the maximum and minimum value of kappa.
            # Used to draw the graph with less k lines
            temp_kappa_list = []
            for i in range(len(kappa_list)):
                temp_kappa_list.append(kappa_list[i] + std_kappa_list[i])
            self.max_kappa = max(temp_kappa_list)
            temp_kappa_list = []
            for i in range(len(kappa_list)):
                temp_kappa_list.append(kappa_list[i] - std_kappa_list[i])
            self.min_kappa = min(temp_kappa_list)
            # get the limit of the k lines close to the kappa points
            i = 2
            kappa = 1
            step = 0.1
            kappa_start_pos = 2
            kappa_end_pos = len(header)
            while True:
                if self.max_kappa > kappa:
                    kappa_start_pos = max(2, i - 3)
                    break
                i += 1
                kappa -= step
                if kappa == step:
                    step /= 10
                if i >= len(header):
                    kappa_start_pos = len(header)
                    break
            i = 2
            kappa = 1
            step = 0.1
            while True:
                if self.min_kappa > kappa:
                    kappa_end_pos = min(i + 3, len(header))
                    break
                i += 1
                kappa -= step
                if kappa == step:
                    step /= 10
                if i >= len(header):
                    kappa_end_pos = len(header)
                    break
        if self.kappa_ax is None:
            figure, ax = plt.subplots(facecolor=settings.GRAPH_BACKGROUND_COLOR)
            ax.axes.set_frame_on(False)
            ax.grid(True, which='both', color=GRID_COLOR)
            ax.axhline(0, color=AX_LINE_COLOR, linewidth=4)
            ax.axvline(0, color=AX_LINE_COLOR, linewidth=4)
            ax.tick_params(color=AX_TICK_COLOR, which='both', labelcolor=AX_TICK_COLOR, labelsize=AX_TICK_SIZE)
            ax.axes.set_ylim([0.1, 1.5])
            ax.axes.set_xlim([10, 200])
            ax.set_xlabel("Dry diameter(nm)", color=LABEL_COLOR, size=LABEL_SIZE)
            ax.set_ylabel("Super Saturation(%)", color=LABEL_COLOR, size=LABEL_SIZE)
            if (self.is_show_all_k_points):
                ax.set_title("Activation Diameter for all Kappa points and Lines of Constant Kappa (K)",
                             color=TITLE_COLOR, size=TITLE_SIZE)
            else:
                ax.set_title("Activation Diameter for average Kappa Points and Lines of Constant Kappa (K)",
                             color=TITLE_COLOR, size=TITLE_SIZE)
            figure.canvas.mpl_connect('pick_event', self.on_pick_kappa_points)
            # graph the k lines
            for i in range(kappa_start_pos, kappa_end_pos):
                y = self.klines_data[header[i]]
                ax.loglog(diameter_list, y, label=str(header[i]), linewidth=4)
            # Graph all the kappa points
            self.kappa_points = ax.scatter(k_points_x_list, k_points_y_list, s=300, c=all_k_points_color_list, picker=5,
                                           label="kappa points")
            x = k_points_x_list[self.current_kappa_point_index]
            y = k_points_y_list[self.current_kappa_point_index]
            # graph the current selection point
            self.current_kappa_point, = ax.plot(x, y, 'o', color=KAPPA_CURRENT_SELECTION_COLOR,
                                                mew=0.5, ms=18, label="current selection")
            handles, labels = ax.get_legend_handles_labels()
            legend = ax.legend(handles, labels)
            self.kappa_figure = figure
            self.kappa_ax = ax

        else:
            if (self.is_show_all_k_points):
                self.kappa_ax.set_title("Activation Diameter for all Kappa points and Lines of Constant Kappa (K)",
                                        color=TITLE_COLOR, size=TITLE_SIZE)
            else:
                self.kappa_ax.set_title("Activation Diameter for average Kappa Points and Lines of Constant Kappa (K)",
                                        color=TITLE_COLOR, size=TITLE_SIZE)
            self.kappa_points.set_offsets(numpy.c_[k_points_x_list, k_points_y_list])
            self.kappa_points.set_color(all_k_points_color_list)
            x = k_points_x_list[self.current_kappa_point_index]
            y = k_points_y_list[self.current_kappa_point_index]
            self.current_kappa_point.set_xdata(x)
            self.current_kappa_point.set_ydata(y)

    def on_key_release_kappa_graph(self, event):
        """
        When keys are pressed. Used to select current kappa points
        :param event: key press event
        :return:
        """
        key = event.key()
        if self.current_kappa_point_index is None:
            self.current_kappa_point_index = len(self.kappa_points_data_list) - 1
        # left arrow key
        elif (key == 16777234):
            self.current_kappa_point_index = min(len(self.kappa_points_data_list) - 1,
                                                 self.current_kappa_point_index + 1)
        # right arrow key
        elif (key == 16777236):
            self.current_kappa_point_index = max(0, self.current_kappa_point_index - 1)
        # up arrow key
        elif (key == 16777235):
            curr_ss = self.kappa_points_data_list[self.current_kappa_point_index][1]
            for i in range(len(self.kappa_points_data_list) - 1, -1, -1):
                if self.kappa_points_data_list[i][1] > curr_ss:
                    self.current_kappa_point_index = i
                    break
        # down arrow key
        elif (key == 16777237):
            curr_ss = self.kappa_points_data_list[self.current_kappa_point_index][1]
            for i in range(len(self.kappa_points_data_list)):
                if self.kappa_points_data_list[i][1] < curr_ss:
                    self.current_kappa_point_index = i
                    break
        self.update_kappa_info_and_graph()

    def on_pick_kappa_points(self, event):
        """
        When mouse is clicked on a kappa point, change to focus that kappa point
        :param event:
        :return:
        """
        self.current_kappa_point_index = event.ind[0]
        self.update_kappa_info_and_graph()

    def on_select_kappa_points_through_info_table(self, row):
        """
        When a kappa point is selected by clicking on the info section representing the point
        :param row:
        :return:
        """
        self.current_kappa_point_index = row - 1
        self.update_kappa_info_and_graph()

    def toggle_exclude_include_kappa_point(self):
        """
        Toggle between including and excluding the current kappa point. If it is already excluded, then include it.
        Otherwise, exclude it.
        :return:
        """
        ss = self.kappa_points_data_list[self.current_kappa_point_index][1]
        dp = self.kappa_points_data_list[self.current_kappa_point_index][0]
        if self.kappa_points_is_included_list[(dp, ss)]:
            self.kappa_points_is_included_list[(dp, ss)] = False
        else:
            self.kappa_points_is_included_list[(dp, ss)] = True
        self.update_kappa_info_and_graph()

    def export_to_csv(self):
        """
        Export the resulting data to csv.
        :return:
        """
        file_name = "kappa_" + self.experiment_date.replace("/", ".") + ".xlsx"
        kappa_dict = self.kappa_calculate_dict
        key_list = self.kappa_points_data_list
        usability_list = self.kappa_points_is_included_list
        data = []
        for a_key in key_list:
            a_row = []
            if kappa_dict[a_key[1]]:
                for a_scan_data in kappa_dict[a_key[1]]:
                    if a_scan_data[0] == a_key[0]:
                        a_row.append(a_key[1])
                        a_row.append(a_scan_data[0])
                        a_row.append(a_scan_data[1])
                        a_row.append(a_scan_data[2])
                        a_row.append(a_scan_data[3])
                        break
                data.append(a_row)
        data = numpy.array(data)
        df = pandas.DataFrame(data, columns=["Super Saturation(%)", "dp(nm)", "K/app", "K/ana", "deviation(%"])
        df.to_excel(file_name, sheet_name='Kappa', index=False)
        self.view.show_error_dialog("Export to " + file_name + " is successful!")

    ##############################################
    #
    # COMMON USES
    #
    ##############################################

    def prepare_scan_data(self):
        """
        Prepare the data of each scan for various processing and graphing steps.
        :return:
        """
        try:
            normalized_concentration_list = []
            startPoint = self.current_scan * self.scan_duration
            endPoint = (self.current_scan + 1) * self.scan_duration
            self.particle_diameter_list = list(self.processed_data[startPoint:endPoint]['dp'])
            self.ccn_list = list(self.processed_data[startPoint:endPoint]['CCNC Count'])
            self.cn_list = list(self.processed_data[startPoint:endPoint]['SMPS Count'])
            self.temp1 = list(self.processed_data[startPoint:endPoint]['T1'])
            self.temp2 = list(self.processed_data[startPoint:endPoint]['T2'])
            self.temp3 = list(self.processed_data[startPoint:endPoint]['T3'])
            self.cn_list = [x * 1 / (self.flow_rate * 1000 / 60) for x in self.cn_list]
            self.diameter_midpoint_list = []
            self.ccn_cn_ratio_list = []
            for i in range(len(self.ccn_list)):
                self.ccn_cn_ratio_list.append(self.ccn_list[i] / self.cn_list[i])
            self.drop_size_list = list(self.processed_data[startPoint:endPoint]['Ave Size'])
            for i in range(2, len(self.normalized_concentration_list)):
                self.diameter_midpoint_list.append(self.normalized_concentration_list[i][0])
                normalized_concentration_list.append(self.normalized_concentration_list[i][self.current_scan + 1])
            self.ccn_normalized_list = normalize_list(normalized_concentration_list)
            self.super_saturation_rate = self.super_saturation_list[self.current_scan]
            if self.finish_scan_alignment_and_auto_sig_fit:
                self.ccnc_sig_list = self.ccnc_sig_list_list[self.current_scan]
                self.b = self.b_list[self.current_scan]
                self.d = self.d_list[self.current_scan]
                self.c = self.c_list[self.current_scan]
                self.min_dp_asym = self.min_dp_asym_list[self.current_scan]
                self.min_dp = self.min_dp_list[self.current_scan]
                self.max_dp_asym = self.max_dp_asym_list[self.current_scan]
                self.dp50 = self.dp50_list[self.current_scan][0]
                self.dp50_wet = self.dp50_wet_list[self.current_scan]
                self.dp50_plus_20 = self.dp50_plus_20_list[self.current_scan]
                self.dp50_less_count = 0
                self.dp50_more_count = 0
                for i in self.particle_diameter_list:
                    if i >= self.dp50:
                        self.dp50_more_count += 1
                    else:
                        self.dp50_less_count += 1
                self.ccn_cn_sim_list = [0]
                for i in range(1, len(self.particle_diameter_list)):
                    if self.min_dp < self.particle_diameter_list[i] < self.max_dp:
                        n = self.b / (1 + (self.particle_diameter_list[i] / self.d) ** self.c)
                        self.ccn_cn_sim_list.append(n)
                    else:
                        self.ccn_cn_sim_list.append(self.ccn_cn_sim_list[i - 1])
        except:
            self.min_pos_CCNC_list[self.current_scan] = None
            self.min_pos_SMPS_list[self.current_scan] = None
            self.is_usable_for_sigmoid_fit_list[self.current_scan] = False
            self.is_usable_for_kappa_cal_list[self.current_scan] = False

    def shift_data_by_one_second(self, forward=True):
        """
        Shift the ccnc data by one second. Used to manually align the smps and ccnc data.
        :param forward: true when + 1 second, false otherwise
        :return:
        """
        try:
            # can't shift if already pass the alignment step
            if self.finish_scan_alignment_and_auto_sig_fit:
                self.view.show_error_dialog("Can't shift current run. Already passed the shifting phase!")
                return
            # invalid scan. Can't shift
            if not self.is_usable_for_sigmoid_fit_list[self.current_scan] or self.min_pos_SMPS_list[self.current_scan] \
                    is None or self.min_pos_CCNC_list[self.current_scan] is None:
                return
            shift_factor = self.shift_factor_list[self.current_scan]
            if forward is True:
                shift_factor -= 1
            else:
                shift_factor += 1

            self.shift_factor_list[self.current_scan] = shift_factor
            start_time = self.scan_position_in_data[self.current_scan] + shift_factor
            end_time = start_time + self.scan_duration
            new_ccn = list(self.raw_data.iloc[start_time:end_time, 4])
            new_ave_size = list(self.raw_data.iloc[start_time:end_time, 5])
            new_t1 = list(self.raw_data.iloc[start_time:end_time, 6])
            new_t2 = list(self.raw_data.iloc[start_time:end_time, 7])
            new_t3 = list(self.raw_data.iloc[start_time:end_time, 8])
            data_start_time = self.scan_duration * self.current_scan
            data_end_time = data_start_time + self.scan_duration
            for i in range(self.scan_duration):
                self.processed_data.iat[data_start_time + i, 4] = new_ccn[i]
                self.processed_data.iat[data_start_time + i, 5] = new_ave_size[i]
                self.processed_data.iat[data_start_time + i, 6] = new_t1[i]
                self.processed_data.iat[data_start_time + i, 7] = new_t2[i]
                self.processed_data.iat[data_start_time + i, 8] = new_t3[i]
            self.prepare_scan_data()
            self.draw_concentration_over_scan_time_graph()
            self.draw_ccn_cn_ratio_over_diameter_graph()
            self.draw_temperature_graph()
            self.view.update_scan_information()
        except:
            self.view.show_error_dialog("Can't shift current run. You should disable this run!")

    def change_scan_status(self):
        """
        change the status of a scan from disabled to enabled. Only work with cetain types of scan
        :return:
        """
        # if the scan is marked as untouched, then clear that mark
        for i in range(len(self.unfinished_sigmoid_fit_scans_list)):
            if self.unfinished_sigmoid_fit_scans_list[i] == self.current_scan:
                self.unfinished_sigmoid_fit_scans_list = self.unfinished_sigmoid_fit_scans_list[:i] \
                                                         + self.unfinished_sigmoid_fit_scans_list[i + 1:]
                self.update_view()
                return
        # is aligning the CCNC and SMPS data of each scan
        if self.finish_scan_alignment_and_auto_sig_fit is False:
            if self.is_usable_for_sigmoid_fit_list[self.current_scan] is True:
                self.is_usable_for_sigmoid_fit_list[self.current_scan] = False
            elif self.min_pos_CCNC_list[self.current_scan] and self.min_pos_CCNC_list[self.current_scan]:
                self.is_usable_for_sigmoid_fit_list[self.current_scan] = True
            else:
                self.view.show_error_dialog("You can't enable this scan! This scan is not usable!")
            self.update_view()
        # is checking sigmoid fit
        else:
            if self.is_usable_for_sigmoid_fit_list[self.current_scan]:
                if self.is_usable_for_kappa_cal_list[self.current_scan]:
                    self.is_usable_for_kappa_cal_list[self.current_scan] = False
                else:
                    self.is_usable_for_kappa_cal_list[self.current_scan] = True
            else:
                self.view.show_error_dialog("You can't enable this scan! This scan is not usable!")
            self.update_view()

    def switch_to_scan(self, scan):
        """
        Switch to a new peak
        :param scan: the new peak
        """
        if scan != self.current_scan:
            self.current_scan = scan
            self.update_view()

    def update_view(self):
        """
        Update the content of the view after processing
        :return:
        """
        self.prepare_scan_data()
        self.draw_concentration_over_scan_time_graph()
        if self.finish_scan_alignment_and_auto_sig_fit:
            self.draw_all_scans_alignment_summary_graph()
            self.draw_complete_sigmoid_graph()
            self.view.update_scan_information_after_sigmoid_fit()
        else:
            self.draw_ccn_cn_ratio_over_diameter_graph()
            self.draw_temperature_graph()
            self.view.update_scan_information()
        self.view.centralWidget().setFocus()

    def reset(self):
        """
        Reset the entire UI to work with a new dataset
        :return:
        """
        gc.collect()
        # clear memory of concentration over scan time graph
        if self.concentration_over_scan_time_figure:
            self.concentration_over_scan_time_figure.clf()
            plt.close(self.concentration_over_scan_time_figure)
            self.concentration_over_scan_time_figure = None
        # clear memory of ccn/cn graph
        if self.ccn_cn_ratio_figure:
            self.ccn_cn_ratio_figure.clf()
            plt.close(self.ccn_cn_ratio_figure)
            self.ccn_cn_ratio_figure = None
        # clear memory of temperature graph
        if self.temperature_figure:
            self.temperature_figure.clf()
            plt.close(self.temperature_figure)
            self.temperature_figure = None
        # clear memory of sigmoid graph
        if self.sigmoid_fit_figure:
            self.sigmoid_fit_figure.clf()
            plt.close(self.sigmoid_fit_figure)
            self.sigmoid_fit_figure = None
        # clear memory of scan alignment summary grpah
        if self.all_scans_alignment_figure:
            self.all_scans_alignment_figure.clf()
            plt.close(self.all_scans_alignment_figure)
            self.all_scans_alignment_figure = None
        # clear memory of kappa figure
        if self.kappa_figure:
            self.kappa_figure.clf()
            plt.close(self.kappa_figure)
            self.kappa_figure = None
        self.__init__(self.view)
        gc.collect()

    # ---------------------progress bar-------------------------

    def cancel_progress_bar(self):
        self.cancelling_progress_bar = True

    def start_or_restart(self, data_files):
        # reset the view
        # self.view.reset()
        # reset the controller
        # self.reset()
        # take in new data and rock on!
        self.data_files = data_files
        # we can't do anything without having the Counts2ConcConv constant
        self.counts_to_conc_conv = self.view.get_counts_to_conc_conv()
        # great, now that we have the files, let's parse them
        started_fn = self.view.init_progress_bar
        finish_fn = self.view.close_progress_bar
        progress_fn = self.view.update_progress_bar
        def parse_new_data_thread_sequence(progress_start, progress_update):
            progress_start.emit("Reading in SMPS and CCNC data files")
            self.parse_files()
            progress_update.emit(5)
            self.get_scan_timestamps()
            progress_update.emit(15)
            self.get_normalized_concentration()
            progress_update.emit(35)
            self.get_smps_counts()
            progress_update.emit(65)
            self.get_ccnc_counts()
            progress_update.emit(85)
            self.do_basic_trans()
            progress_update.emit(95)
            self.pre_align_sanity_check()
            progress_update.emit(100)
        worker=run_new_thread(parse_new_data_thread_sequence, started_fn,finish_fn,progress_fn)
        self.thread_pool.start(worker)
        #
        while len(self.scans) == 0:
            time.sleep(0.1)
        self.view.show_alignment_dialog()

    def set_smooth_method(self,method_index):
        if 0 <= method_index < len(CONST.smooth_algos):
            self.smooth_method = CONST.smooth_algos[method_index]
        else:
            self.smooth_method = CONST.smooth_algos[0]

    def set_base_shift_factor(self,n):
        self.base_shift_factor =n

def main():
    example_files = [r"C:\Users\KhaiNguyen\OneDrive\chemics_on_work\Examples\smps.txt",
                     r"C:\Users\KhaiNguyen\OneDrive\chemics_on_work\Examples\ccnc.csv"]
    data_files = os.listdir(r"C:\Users\KhaiNguyen\OneDrive\chemics_on_work\data")
    # os.chdir("../../data")
    # for i in range(len(data_files)):
    #     data_files[i] = os.path.join(os.getcwd(), data_files[i])
    # controller.files = example_files
    # controller.parse_files()
    # controller.get_scan_timestamps()
    # controller.get_normalized_concentration()
    # controller.get_smps_counts()
    # controller.get_ccnc_counts()
    # controller.align_smps_ccnc_data()
    print()
    # view.show_ui()

#
# if __name__ == '__main__':
#     main()
