from Scan import *
import cPickle

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
        self.cancelling_progress_bar = False
        # smps and ccnc files
        self.data_files = None
        # ccnc = Cloud Condensation Nuclei Counter
        self.ccnc_data = None
        # smps = Scanning Mobility Particle Sizer
        self.smps_data = None
        # date of experiment
        self.experiment_date = None
        # smoothing method
        self.smooth_method = CONST.smooth_algos[0]
        # base shift factor. Very useful for auto alignment
        self.base_shift_factor = 0
        # current scan
        self.curr_scan_index = 0
        # limits of the variable b
        self.b_limits = [0.5, 1.5]
        # threshold for top assym
        self.asym_limits = [0.75, 1.5]
        # the current step
        self.stage = "init"
        # other trivial things
        self.save_name = None

        # variables for calculating kappa.
        self.sigma = 0.072
        self.temp = 298.15
        self.dd_1 = 280
        self.i_kappa_1 = 0.00567
        self.dd_2 = 100
        self.i_kappa_2 = 0.6
        self.solubility = 0.03
        self.kappa_excel = None
        self.kappa_calculate_dict = {}
        self.alpha_pinene_dict = {}
        # whether a kappa point is included in calculating the average k
        self.is_valid_kappa_points = {}  # format of key is (dp,ss)

    def init_attributes(self):
        self.scans = []
        self.counts_to_conc_conv = 0
        self.data_files = None
        self.ccnc_data = None
        self.smps_data = None
        self.experiment_date = None
        self.smooth_method = CONST.smooth_algos[0]
        self.base_shift_factor = 0
        self.curr_scan_index = 0
        self.b_limits = [0.5, 1.5]
        self.asym_limits = [0.75, 1.5]
        self.stage = "init"
        self.kappa_calculate_dict = {}
        self.alpha_pinene_dict = {}
        # whether a kappa point is included in calculating the average k
        self.is_valid_kappa_points = {}  # format of key is (dp,ss)
        self.save_name = None

    ##############################################
    #
    # PARSE FILES AND MATCH SMPS AND CCNC DATA
    #
    ##############################################

    def load_project(self, project_file):
        with open(project_file, 'rb') as handle:
            (self.scans, self.counts_to_conc_conv, self.data_files, self.ccnc_data, self.smps_data,
             self.experiment_date, self.smooth_method, self.base_shift_factor, self.b_limits,
             self.asym_limits, self.kappa_calculate_dict, self.alpha_pinene_dict, self.stage,
             self.is_valid_kappa_points,
             self.save_name) = cPickle.load(handle)
        # got to reset the view first
        self.view.reset_view()
        # if got to kappa, go straight to the kappa view
        if self.stage == "kappa":
            self.view.switch_to_kappa_view()
            return
        if self.stage == "sigmoid":
            self.view.show_sigmoid_docker()
        self.view.update_experiment_info()
        self.switch_to_scan(0)

    def save_project(self):
        if self.save_name is None:
            self.view.save_project_as()
        else:
            to_save = (self.scans, self.counts_to_conc_conv, self.data_files, self.ccnc_data, self.smps_data,
                       self.experiment_date, self.smooth_method, self.base_shift_factor, self.b_limits,
                       self.asym_limits, self.kappa_calculate_dict,  self.alpha_pinene_dict , self.stage,
                       self.is_valid_kappa_points,
                       self.save_name)
            with open(self.save_name, 'wb') as handle:
                cPickle.dump(to_save, handle, protocol=cPickle.HIGHEST_PROTOCOL)
                print 'saved'

    def parse_files(self):
        smps_txt_files = []
        ccnc_csv_files = []
        # acquire the smps and ccnc files from the input files
        for a_file in self.data_files:
            if a_file.lower().endswith('.txt'):
                smps_txt_files.append(a_file)
            elif a_file.lower().endswith('.csv'):
                ccnc_csv_files.append(a_file)
        smps_txt_files = [str(x) for x in smps_txt_files]
        smps_txt_files = smps_txt_files[0]
        ccnc_csv_files = [str(x) for x in ccnc_csv_files]
        self.ccnc_data = process_csv_files(ccnc_csv_files)
        self.smps_data = process_text_files(smps_txt_files)
        self.experiment_date = self.ccnc_data[0]
        self.ccnc_data = self.ccnc_data[1]

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
            a_scan = Scan()
            self.scans.append(a_scan)
            start_time = datetime.strptime(scan_start_times[i], "%H:%M:%S")
            end_time = start_time + timedelta(seconds=scan_duration)
            a_scan.set_start_time(start_time)
            a_scan.set_end_time(end_time)
            a_scan.set_up_time(scan_up_time)
            a_scan.set_down_time(scan_down_time)
            a_scan.set_duration(scan_duration)
            a_scan.set_counts_2_conc(self.counts_to_conc_conv)

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
            for j in range(len(self.scans)):
                a_scan = self.scans[j]
                a_scan.add_to_diameter_midpoints(self.smps_data[i][0])
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
            if are_floats_equal(curr_time, target_time) or curr_line_index == end_line_index:
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
        # todo: outliers detection
        # todo: fix all wrong shift factors
        pass

    def preview_scans(self, timer):
        # becase there is an inherent delay time when switch
        timer = max(0, timer - 0.5)
        temp_curr_scan_index = self.curr_scan_index
        self.view.init_progress_bar("Previewing all scans...")
        for i in range(len(self.scans)):
            self.view.update_progress_bar(100 * i // len(self.scans))
            self.switch_to_scan(i)
            if self.view.progress_dialog.wasCanceled():
                break
            time.sleep(timer)
        self.switch_to_scan(temp_curr_scan_index)
        self.view.close_progress_bar()

    def reset_project(self):
        self.start(self.data_files)

    def start(self, data_files):
        # reset the attributes
        self.init_attributes()
        # got to reset the view first
        self.view.reset_view()
        # take in new data and rock on!
        self.data_files = data_files
        # we can't do anything without having the Counts2ConcConv constant
        self.counts_to_conc_conv = self.view.get_counts_to_conc_conv()
        if self.counts_to_conc_conv is None:
            return
        # great, now that we have the files, let's parse them
        self.view.init_progress_bar("Reading in SMPS and CCNC data files")
        self.parse_files()
        self.view.update_progress_bar(5)
        self.get_scan_timestamps()
        self.view.update_progress_bar(15)
        self.get_normalized_concentration()
        self.view.update_progress_bar(35)
        self.get_smps_counts()
        self.view.update_progress_bar(65)
        self.get_ccnc_counts()
        self.view.update_progress_bar(85)
        self.do_basic_trans()
        self.view.update_progress_bar(95)
        self.pre_align_sanity_check()
        self.view.update_progress_bar(100)
        self.view.close_progress_bar()
        # update some basic info first
        self.view.update_experiment_info()
        self.view.show_auto_or_manual_question_dialog()
        self.stage = "align"
        # update the menu
        self.view.set_menu_bar_by_stage()

    def align_smps_ccnc_data(self):
        shift_factor = self.base_shift_factor
        self.view.init_progress_bar("Aligning SMPS and CCNC data...")
        for i in range(len(self.scans)):
            self.view.update_progress_bar(100 * (i + 1) // len(self.scans))
            shift_factor = self.scans[i].align_smps_ccnc_data(shift_factor)
            self.scans[i].generate_processed_data()
        self.view.close_progress_bar()
        self.post_align_sanity_check()
        self.switch_to_scan(0)

    def prepare_data_for_manual_inputs(self):
        self.view.init_progress_bar("Processing SMPS and CCNC data...")
        for i in range(len(self.scans)):
            self.view.update_progress_bar(100 * (i + 1) // len(self.scans))
            self.scans[i].set_shift_factor(0)
            self.scans[i].generate_processed_data()
        self.view.close_progress_bar()
        self.switch_to_scan(0)

    ##############################################
    #
    # CORRECT CHARGES AND FIT SIGMOID
    #
    ##############################################
    def correct_charges(self):
        # Correcting charges
        self.view.init_progress_bar("Correcting charges...")
        for i in range(len(self.scans)):
            # for i in range(len(self.scans)):
            self.view.update_progress_bar(100 * (i + 1) // len(self.scans))
            self.scans[i].correct_charges()
        self.view.close_progress_bar()
        self.view.show_sigmoid_docker()
        self.stage = "sigmoid"
        self.switch_to_scan(0)

    def auto_fit_sigmoid(self):
        self.view.init_progress_bar("Auto fitting sigmoid...")
        for i in range(len(self.scans)):
            self.view.update_progress_bar(100 * (i + 1) // len(self.scans))
            self.scans[i].cal_params_for_sigmoid_fit()
            self.scans[i].fit_sigmoids()
        self.switch_to_scan(0)

    ##############################################
    #
    # KAPPA CALCULATION
    #
    ##############################################
    def cal_kappa(self):
        self.calculate_all_kappa_values()
        self.calculate_average_kappa_values()
        self.stage = "kappa"
        self.view.switch_to_kappa_view()

    def calculate_all_kappa_values(self):
        """
        Calculate the kappa values for every scan.
        """
        if self.kappa_excel is None:
            self.kappa_excel = pandas.read_csv("kCal.csv", header=None)
        lookup = self.kappa_excel
        a_param = 0.00000869251 * self.sigma / self.temp
        # asc = (exp(sqrt(4 * a_param ** 3 / (27 * self.i_kappa_1 * (self.dd_1 * 0.000000001) ** 3))) - 1) * 100
        # Calculate each kappa
        ss_and_dps = []
        for i in range(len(self.scans)):
            scan = self.scans[i]
            if not scan.is_valid() or scan.true_super_sat is None:
                continue
            ss = scan.true_super_sat
            for j in range(len(scan.dps)):
                dp_50 = scan.dps[j][0]
                ss_and_dps.append([ss, dp_50])
        for i in range(len(ss_and_dps)):
            ss = float(ss_and_dps[i][0])
            dp_50 = float(ss_and_dps[i][1])
            row_index = int(math.floor(dp_50 - 9))
            match_row = list(lookup.iloc[row_index][1:])
            value_row = list(lookup.iloc[0][1:])
            a = get_correct_num(match_row, ss)
            c_index = a[1]
            a = a[0]
            if c_index != (len(match_row) - 1):
                b = match_row[c_index + 1]
                c = value_row[c_index]
                d = value_row[c_index + 1]
            else:
                c = value_row[c_index]
                b = 0
                d = 0
            apparent_kappa = (ss - (a - (a - b) / (c - d) * c)) / ((a - b) / (c - d))
            analytic_kappa = (4 * a_param ** 3) / (27 * (dp_50 * 0.000000001) ** 3 * log(ss / 100 + 1) ** 2)
            deviation_percentage = (apparent_kappa - analytic_kappa) / apparent_kappa * 100
            if ss in self.kappa_calculate_dict.keys():
                self.kappa_calculate_dict[ss].append([dp_50, apparent_kappa, analytic_kappa, deviation_percentage])
            else:
                self.kappa_calculate_dict[ss] = ([[dp_50, apparent_kappa, analytic_kappa, deviation_percentage]])
            self.is_valid_kappa_points[(dp_50, ss)] = True

    def calculate_average_kappa_values(self):
        """
        Calculate the kappa values for each super saturation percentage. The values are average of all scans with the
        same super saturation
        :return:
        """
        self.alpha_pinene_dict = {}
        for a_key in self.kappa_calculate_dict.keys():
            a_scan = self.kappa_calculate_dict[a_key]
            temp_dp50_list = []
            dp_50s = []
            apparent_kappas = []
            analytical_kappas = []
            mean_of_stds = []
            for aSS in a_scan:
                dp_50s.append(aSS[0])
                if self.is_valid_kappa_points[(aSS[0], a_key)]:
                    temp_dp50_list.append(aSS[0])
                    apparent_kappas.append(aSS[1])
                    analytical_kappas.append(aSS[2])
                    mean_of_stds.append(aSS[3])
            mean_dp = average(temp_dp50_list)
            std_dp = numpy.std(temp_dp50_list)
            mean_app = average(apparent_kappas)
            std_app = numpy.std(apparent_kappas)
            mean_ana = average(analytical_kappas)
            std_ana = numpy.std(analytical_kappas)
            mean_dev = average(mean_of_stds)
            dev_mean = (mean_app - mean_ana) / mean_app * 100
            self.alpha_pinene_dict[a_key] = (
                mean_dp, std_dp, mean_app, std_app, mean_ana, std_ana, mean_dev, dev_mean, dp_50s)
            print()


    def export_project_data(self):
        """
        Export the resulting data to csv.
        :return:
        """
        file_name = "kappa_" + self.experiment_date.replace("/", ".") + ".csv"
        data = []
        for a_key in self.kappa_calculate_dict.keys():
            a_scan = self.kappa_calculate_dict[a_key]
            for aSS in a_scan:
                if self.is_valid_kappa_points[(aSS[0], a_key)]:
                    a_row = [a_key] + aSS + ["Included point"]
                else:
                    a_row = [a_key] + aSS + ["Excluded point"]
                data.append(a_row)
        df = pandas.DataFrame(numpy.asarray(data), columns=["Super Saturation(%)", "dp(nm)", "K/app", "K/ana",
                                                            "deviation(%","Status"])
        df.to_csv(file_name, index=False)
        self.view.show_information_message(title = "Export Data", text= "Export to " + file_name + " successful!")

    ##############################################
    #
    # COMMON USES
    #
    ##############################################

    def switch_to_scan(self, index):
        self.curr_scan_index = index
        self.view.update_scan_info_and_graphs()
        self.view.set_menu_bar_by_stage()

    def set_smooth_method(self, method_index):
        if 0 <= method_index < len(CONST.smooth_algos):
            self.smooth_method = CONST.smooth_algos[method_index]
        else:
            self.smooth_method = CONST.smooth_algos[0]

    def set_base_shift_factor(self, n):
        self.base_shift_factor = n

    def set_scan_index(self, n):
        self.curr_scan_index = n

    def set_save_name(self, name):
        self.save_name = name

    def set_sigma(self, sigma):
        self.sigma = sigma

    def set_temp(self, temp):
        self.temp = temp

    def set_dd_1(self, dd_1):
        self.dd_1 = dd_1

    def set_dd_2(self, dd_2):
        self.dd_2 = dd_2

    def set_i_kappa_1(self, value):
        self.i_kappa_1 = value

    def set_i_kapp_2(self, value):
        self.i_kappa_2 = value

    def set_solubility(self, value):
        self.solubility = value

    def set_kappa_point_state(self, ss, dp, state):
        self.is_valid_kappa_points[(dp, ss)] = state
        self.calculate_average_kappa_values()
        self.view.update_kappa_graph()
