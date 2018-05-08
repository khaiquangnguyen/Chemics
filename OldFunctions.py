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
                a_param / aNum))

    # Calculate the third colum
    dList2 = [self.dd_1 * 0.000000001]
    for i in range(1000):
        dList2.append(dList2[-1] * 1.005)

    # Calculate the fourth column
    sList2 = []
    firstNum = dList2[0]
    for i in range(len(dList1)):
        aNum = dList2[i]
        sList2.append(
            (aNum ** 3 - firstNum ** 3) / (aNum ** 3 - firstNum ** 3 * (1 - self.i_kappa_1)) * math.exp(
                a_param / aNum))

    dList3 = [self.dd_2 * 0.000000001]
    for i in range(1000):
        dList3.append(dList3[-1] * 1.005)

    kappaList = []
    firstNum = dList3[0]
    for i in range(len(dList3)):
        aNum = dList3[i]
        kappaList.append(min((aNum ** 3 / firstNum ** 3 - 1) * self.solubility, 1) * self.i_kappa_2)

    sList3 = []
    firstNum = dList3[0]
    for i in range(len(dList1)):
        aNum = dList3[i]
        if kappaList[i] == 0:
            sList3.append(0)
        else:
            sList3.append(
                (aNum ** 3 - firstNum ** 3) / (aNum ** 3 - firstNum ** 3 * (1 - kappaList[i])) * math.exp(
                    a_param / aNum))
    self.sc = (max(sList2) - 1) * 100
    self.sc2 = (max(sList3) - 1) * 100
    scCalcs = True

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


# def find_ref_index_ccnc(ccnc_list,potential_loc):
#     # we are working with a lot of assumptions here
#     # The first one is that the two ref values must be quite close to each other
#     # find absolute max
#     # the potential location is the sum of the smps local minimum and the base shift factor calculated from other
#     # scans
#     left_max = numpy.argmax(ccnc_list)
#     # find the right list which contains the lower point.
#     right_list = ccnc_list[left_max:]
#     # smooth the hell out of the data
#     right_list = smooth(right_list)
#     #find min
#     potential_mins = scipy.signal.argrelmin(right_list, order=2)
#     # if find no min, return error
#     # plt.plot(right_list)
#     # plt.show()
#     if len(potential_mins[0]) == 0:
#         return None
#     print potential_mins
#     # now, compare among the potential mins to see if which one has the closest value to the smps low
#     for i in range(len(potential_mins[0])):
#         ref_index = left_max + potential_mins[0][i]
#         # if agree with the potential loc
#         if abs(ref_index - potential_loc) <= 2:
#             return ref_index
#     return None

def show_error_dialog(self, error_message='Unknown Error!'):
    """
    Show the error message
    :param errorMessage: The message to show in the error message
    """
    self.progress_dialog.reset_view()
    warning = QMessageBox()
    warning.setIcon(QMessageBox.Warning)
    warning.setText(error_message)
    warning.exec_()


def get_min_index_ccnc_old(a_list,scan_up_time = 0):
    """
    get the position of the smallest value of aParam list
    :param processed_data: the processed_data to process
    :return: the index of the smallest value. -1 if the peak is not usable
    """
    a_list = numpy.asarray(a_list)
    # TODO: improve this method
    first_max = 0
    max_pos = 0
    if sum(a_list) == 0:
        return -1
    max_pos = peakutils.indexes(a_list, thres=0.5, min_dist=len(a_list) / 10)
    if len(max_pos) == 0:
        return -1
    else:
        max_pos = max_pos[0]
    first_max = a_list[max_pos]

    # Get the second maximum
    second_max = 0
    second_max_dis = 0
    second_max_pos = 0
    for i in range(max_pos + 2, len(a_list)):
        max_dis = 0
        if a_list[i] < first_max * 0.1:
            continue
        for j in range(1, i - max_pos):
            if a_list[i] <= a_list[i - j]:
                break
            else:
                max_dis += 1
        if max_dis >= second_max_dis:
            second_max = a_list[i]
            second_max_pos = i
            second_max_dis = max_dis

    # Check if the two peaks are actually usable
    # Get the minimum between two peaks
    min = first_max
    min_pos = max_pos
    for i in range(max_pos, second_max_pos):
        if a_list[i] <= min:
            min = a_list[i]
            min_pos = i
        if a_list[i] - min < 10:
            min = a_list[i]
            min_pos = i
    # If the second peak is too small,then not valid
    if second_max < first_max * 0.1:
        return -1
    # if the minimum is too big, the not valid as well
    if min >= first_max * 0.9:
        return -1
    # If either peak is 0
    if second_max == 0 or first_max == 0:
        return -1
    return min_pos


def create_temperature_graph(self, new_figure=None):
    try:
        if new_figure is None:
            figure = plt.figure(facecolor=settings.graphBackgroundColor)
        else:
            figure = new_figure
            figure.clf()
            plt.figure(figure.number)
        plt.axes(frameon=False)
        plt.grid(color='0.5')
        plt.axhline(0, color='0.6', linewidth=4)
        plt.axvline(0, color='0.6', linewidth=4)
        plt.gca().tick_params(axis='x', color='1', which='both', labelcolor="0.6")
        plt.gca().tick_params(axis='y', color='1', which='both', labelcolor="0.6")
        plt.gca().yaxis.label.set_color('0.6')
        plt.gca().xaxis.label.set_color('0.6')

        x = range(self.scan_duration)
        minY = min(min(self.temp1), min(self.temp2), min(self.temp3)) - 3
        maxY = max(max(self.temp1), max(self.temp2), max(self.temp3)) + 3
        plt.gca().axes.set_ylim([minY, maxY])
        plt.plot(x, self.temp1, linewidth=5, color='#EF5350', label="T1")
        plt.plot(x, self.temp2, linewidth=5, color='#2196F3', label="T2")
        plt.plot(x, self.temp3, linewidth=5, color='#1565C0', label="T3")
        plt.xlabel("Scan time(s)")
        plt.ylabel("Temperature")
        handles, labels = plt.gca().get_legend_handles_labels()
        legend = plt.legend(handles, labels, loc="upper right", bbox_to_anchor=(1.1, 1.1))
        legend.get_frame().set_facecolor('#9E9E9E')
    except:
        figure = plt.figure(facecolor=settings.graphBackgroundColor)
    finally:
        self.temperature_graph_list.append(plt.gcf())


def create_concentration_over_scan_time_graph(self, new_figure=None):
        try:
            if len(self.cn_list) != len(self.ccn_list) or len(self.cn_list) == 0 or len(self.ccn_list) == 0:
                return
            if new_figure is None:
                figure = plt.figure(facecolor=settings.graphBackgroundColor)
            else:
                figure = new_figure
                figure.clf()
                plt.figure(figure.number)
            plt.axes(frameon=False)
            plt.grid(color='0.5')
            plt.axhline(0, color='0.6', linewidth=4)
            plt.axvline(0, color='0.6', linewidth=4)
            plt.gca().tick_params(axis='x', color='1', which='both', labelcolor="0.6")
            plt.gca().tick_params(axis='y', color='1', which='both', labelcolor="0.6")
            plt.gca().yaxis.label.set_color('0.6')
            plt.gca().xaxis.label.set_color('0.6')
            x = range(self.scan_duration)
            plt.plot(x, self.cn_list, linewidth=4, color='#EF5350', label="CN")
            plt.plot(x, self.ccn_list, linewidth=4, color='#2196F3', label="CCN")
            handles, labels = plt.gca().get_legend_handles_labels()
            legend = plt.legend(handles, labels, loc="upper left", bbox_to_anchor=(0, 0.9))
            legend.get_frame().set_facecolor('#9E9E9E')

            plt.xlabel("Scan time(s)")
            plt.ylabel("Concentration(cm3)")
        except:
            figure = plt.figure(facecolor=settings.graphBackgroundColor)
        finally:
            self.adjusted_graph_list.append(plt.gcf())


def create_ccn_cn_ratio_over_diameter_graph(self, new_figure=None):
        try:
            if new_figure is None:
                figure = plt.figure(facecolor=settings.graphBackgroundColor)
            else:
                figure = new_figure
                figure.clf()
                plt.figure(figure.number)
            plt.axes(frameon=False)
            plt.grid(color='0.5')
            plt.axhline(0, color='0.6', linewidth=2)
            plt.axvline(0, color='0.6', linewidth=4)
            plt.axhline(1, color='0.7', linewidth=2, linestyle='dashed')
            plt.gca().tick_params(axis='x', color='1', which='both', labelcolor="0.6")
            plt.gca().tick_params(axis='y', color='1', which='both', labelcolor="0.6")
            plt.gca().yaxis.label.set_color('0.6')
            plt.gca().xaxis.label.set_color('0.6')
            plt.plot(self.particle_diameter_list, self.ccn_cn_ratio_list, 'o', color="#2196F3", mew=0.5,
                     mec="#0D47A1",
                     ms=9, label="CCN/CN")
            yLim = min(2, max(self.ccn_cn_ratio_list)) + 0.2
            plt.gca().axes.set_ylim([-0.1, yLim])
            handles, labels = plt.gca().get_legend_handles_labels()
            legend = plt.legend(handles, labels, loc="upper right", bbox_to_anchor=(1, 0.9))
            legend.get_frame().set_facecolor('#9E9E9E')
            plt.xlabel("Diameter (nm)")
            plt.ylabel("CCN/CN")
        except:
            figure = plt.figure(facecolor=settings.graphBackgroundColor)
        finally:
            self.dry_diameter_graph_list.append(plt.gcf())


def draw_complete_sigmoid_graph(self, new_figure=None):
        """
        Make complete graph of the dry diameter after optimization and sigmodal fit
        """
        try:
            if not self.usable_for_sigmoid_fit_list[self.current_scan]:
                self.create_ccn_cn_ratio_over_diameter_graph(new_figure)
                return
            if new_figure is None:
                figure = plt.figure(facecolor=settings.graphBackgroundColor)
            else:
                figure = new_figure
                figure.clf()
                plt.figure(figure.number)
            plt.axes(frameon=False)
            plt.grid(color='0.5')
            plt.axhline(0, color='0.6', linewidth=4)
            plt.axvline(0, color='0.6', linewidth=4)
            plt.axhline(1, color='0.7', linewidth=2, linestyle="--")
            plt.gca().tick_params(axis='x', color='1', which='both', labelcolor="0.6")
            plt.gca().tick_params(axis='y', color='1', which='both', labelcolor="0.6")
            plt.gca().yaxis.label.set_color('0.6')
            plt.gca().xaxis.label.set_color('0.6')
            yLim = min(2, max(self.ccnc_sig_list)) + 0.2
            plt.gca().axes.set_ylim([-0.1, yLim])
            plt.plot(self.diameter_midpoint_list, self.ccn_normalized_list, linewidth=4, color='#43A047',
                     label="dN/dLogDp")
            if self.usable_for_kappa_cal_list[self.current_scan] and self.usable_for_sigmoid_fit_list[
                self.current_scan]:
                plt.plot(self.particle_diameter_list, self.ccn_cn_sim_list, linewidth=5, color='#EF5350',
                         label="Sigmodal Fit")
            plt.plot(self.particle_diameter_list, self.ccn_cn_ratio_list, 'o', color="#2196F3", mew=0.5,
                     mec="#1976D2",
                     ms=9, label="CCN/CN")
            plt.plot(self.particle_diameter_list, self.ccnc_sig_list, 'o', color="#1565C0", mew=0.5, mec="#0D47A1",
                     ms=9, label="CCN/CN (Corrected)")
            plt.xlabel("Dry diameter(nm)")
            plt.ylabel("CCN/CN ratio and Normalized dN/dLogDp")
            handles, labels = plt.gca().get_legend_handles_labels()
            legend = plt.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.7, 1.1))
            legend.get_frame().set_facecolor('#9E9E9E')
            self.dry_diameter_graph_list.append(plt.gcf())
        except:
            figure = plt.figure(facecolor=settings.graphBackgroundColor)
            self.dry_diameter_graph_list.append(plt.gcf())


def draw_all_scans_alignment_summary_graph(self):
    """
    A graph of peak alignment, and also allow interaction to select peak to process
    """
    # Prepare the figure
    figure = plt.figure(facecolor=settings.graphBackgroundColor)
    plt.axes(frameon=False)
    plt.grid(color='0.5')
    plt.axhline(0, color='0.6', linewidth=4)
    plt.axvline(0, color='0.6', linewidth=4)
    plt.gca().tick_params(axis='x', color='1', which='both', labelcolor="0.6")
    plt.gca().tick_params(axis='y', color='1', which='both', labelcolor="0.6")
    plt.gca().yaxis.label.set_color('0.6')
    plt.gca().xaxis.label.set_color('0.6')
    figure.canvas.mpl_connect('pick_event', self.on_pick)
    tempSMPSPeakList = []
    tempCCNPeakCList = []
    for i in range(len(self.min_pos_CCNC_list)):
        if self.min_pos_CCNC_list[i] and self.min_pos_SMPS_list[i] and self.usable_for_sigmoid_fit_list[i]:
            tempSMPSPeakList.append(self.min_pos_SMPS_list[i])
            tempCCNPeakCList.append(self.min_pos_CCNC_list[i])
        else:
            # Make up for the null values
            if len(tempSMPSPeakList) > 0:
                tempSMPSPeakList.append(tempSMPSPeakList[-1] + self.scan_duration)
                tempCCNPeakCList.append(tempCCNPeakCList[-1] + self.scan_duration)
            else:
                tempSMPSPeakList.append(10)
                tempCCNPeakCList.append(10)

    x = numpy.asarray(tempSMPSPeakList)
    y = numpy.asarray(tempCCNPeakCList)

    result = scipy.stats.linregress(x, y)
    slope = result[0]
    yIntercept = result[1]

    # Recalculate the position of the smps
    plt.plot(x, x * slope + yIntercept, linewidth=4, color='#43A047', label="Regression line")
    plt.plot(x, y, "o", ms=10, color="#43A047", picker=5, mew=0, label="Minimum")
    slope = ('{0:.4f}'.format(slope))
    yIntercept = ('{0:.4f}'.format(yIntercept))
    textToShow = str(slope) + "* x" + " + " + str(yIntercept)
    self.current_point, = plt.plot(x[0], y[0], 'o', color="#81C784", ms=12, mew=0)
    plt.xlabel("SMPS minumum point")
    plt.ylabel("CCNC minimum point")
    handles, labels = plt.gca().get_legend_handles_labels()
    legend = plt.legend(handles, labels, loc="upper right", bbox_to_anchor=(1, 0.7))
    legend.get_frame().set_facecolor('#9E9E9E')

    if not self.min_compare_graph:
        self.min_compare_graph = plt.gcf()
    else:
        plt.close(self.min_compare_graph)
        self.min_compare_graph = plt.gcf()

    def correct_charges(self):
        try:
            start = timer()
            asymp = 99999
            newList = []
            epsilon = 0.0000000001
            e = scipy.constants.e
            e0 = scipy.constants.epsilon_0
            k = scipy.constants.k
            t = scipy.constants.zero_Celsius + 25
            z = 0.875
            p = 1013
            nair = 0.000001458 * t ** 1.5 / (t + 110.4)
            lambdaAir = 2 * nair / 100 / p / (8 * 28.84 / pi / 8.314 / t) ** 0.5 * 1000 ** 0.5
            coeficientList = [[-0.0003, -0.1014, 0.3073, -0.3372, 0.1023, -0.0105],
                              [-2.3484, 0.6044, 0.48, 0.0013, -0.1553, 0.032],
                              [-44.4756, 79.3772, -62.89, 26.4492, -5.748, 0.5049]]

            # frac0List = calculate_fraction(self.diameterList,0,coeficientList[0])
            frac1List = calculate_fraction(self.particle_diameter_list, 1, coeficientList[1])
            frac2List = calculate_fraction(self.particle_diameter_list, 2, coeficientList[2])
            frac3List = calculate_fraction(self.particle_diameter_list, 3)
            chargeList = []

            for i in self.particle_diameter_list:
                aDList = [0]
                for k in range(1, 4):
                    c = cal_cc(i * 10 ** -9, lambdaAir)
                    dp = 10 ** 9 * FastDpCalculator.find_dp(i * 10 ** -9 / c, lambdaAir, k)
                    aDList.append(dp)
                chargeList.append(aDList)
            # second part of correct charges
            self.cn_fixed_list = self.cn_list[:]
            self.ccn_fixed_list = self.ccn_list[:]
            maxUpperBinBound = (self.particle_diameter_list[-1] + self.particle_diameter_list[-2]) / 2
            lenDpList = len(self.particle_diameter_list)
            for i in range(lenDpList):
                n = lenDpList - i - 1
                moveDoubletCounts = frac2List[n] / (frac1List[n] + frac2List[n] + frac3List[n]) * self.cn_list[n]
                moveTripletCounts = frac3List[n] / (frac1List[n] + frac2List[n] + frac3List[n]) * self.cn_list[n]
                self.cn_fixed_list[n] = self.cn_fixed_list[n] - moveDoubletCounts - moveTripletCounts
                self.ccn_fixed_list[n] = self.ccn_fixed_list[n] - moveDoubletCounts - moveTripletCounts
                if chargeList[n][2] <= maxUpperBinBound:
                    j = lenDpList - 2
                    while (True):
                        upperBinBound = (self.particle_diameter_list[j] + self.particle_diameter_list[j + 1]) / 2
                        lowerBinBound = (self.particle_diameter_list[j] + self.particle_diameter_list[j - 1]) / 2
                        if upperBinBound > chargeList[n][2] >= lowerBinBound:
                            self.cn_fixed_list[j] = self.cn_fixed_list[j] + moveDoubletCounts
                            if chargeList[n][2] < asymp:
                                if self.g_ccn_list[j] > epsilon:
                                    self.ccn_fixed_list[j] = self.ccn_fixed_list[j] + moveDoubletCounts * \
                                                                                      self.g_ccn_list[j] / \
                                                                                      self.g_cn_list[j]
                            else:
                                self.ccn_fixed_list[j] = self.ccn_fixed_list[j] + moveDoubletCounts
                            break
                        j -= 1

                if chargeList[n][3] < maxUpperBinBound:
                    j = lenDpList - 2
                    while (True):
                        upperBinBound = (self.particle_diameter_list[j] + self.particle_diameter_list[j + 1]) / 2
                        lowerBinBound = (self.particle_diameter_list[j] + self.particle_diameter_list[j - 1]) / 2
                        if upperBinBound > chargeList[n][3] >= lowerBinBound:
                            self.cn_fixed_list[j] = self.cn_fixed_list[j] + moveTripletCounts
                            if chargeList[n][3] < asymp:
                                self.ccn_fixed_list[j] = self.ccn_fixed_list[j] + moveTripletCounts * \
                                                                                  self.ccn_list[j] / self.cn_list[j]
                            else:
                                self.ccn_fixed_list[j] = self.ccn_fixed_list[j] + moveTripletCounts
                            break
                        j -= 1
            for i in range(len(self.ccn_fixed_list)):
                if self.ccn_fixed_list[i] / self.cn_fixed_list[i] < -0.01:
                    self.ccn_fixed_list[i] = 0

            self.g_ccn_list = self.ccn_fixed_list[:]
            self.g_cn_list = self.cn_fixed_list[:]
            self.usable_for_kappa_cal_list[self.current_scan] = True

        except:
            raise ScanDataError()

def refitting_sigmoid_line(self, min_dry_diameter, min_dry_diameter_asymptote, max_dry_diameter_asymptote):
    self.move_progress_bar_forward("Refitting sigmoid line to scan #" + str(self.current_scan + 1),
                                   max_value=6)
    self.ccnc_sig_list = []
    try:
        self.prepare_scan_data()
        self.get_parameters_for_sigmoid_fit(min_dry_diameter, min_dry_diameter_asymptote,
                                            max_dry_diameter_asymptote)
        self.fit_sigmoid_line()
        self.move_progress_bar_forward()
        self.min_dp_list[self.current_scan] = self.min_dp
        self.min_dp_asym_list[self.current_scan] = self.min_dp_asym
        self.max_dp_asym_list[self.current_scan] = self.max_dp_asym
        self.b_list[self.current_scan] = self.b
        self.d_list[self.current_scan] = self.d
        self.c_list[self.current_scan] = self.c
        self.dp50_list[self.current_scan] = (self.d, self.super_saturation_list[self.current_scan])
        self.usable_for_kappa_cal_list[self.current_scan] = True
        self.usable_for_sigmoid_fit_list[self.current_scan] = True
        self.update_scan_info_and_graphs()
        self.move_progress_bar_forward(complete=True)
    except:
        self.usable_for_kappa_cal_list[self.current_scan] = False
        self.usable_for_sigmoid_fit_list[self.current_scan] = False
        self.min_dp_list[self.current_scan] = 0
        self.min_dp_asym_list[self.current_scan] = 0
        self.max_dp_asym_list[self.current_scan] = 0
        self.b_list[self.current_scan] = 0
        self.d_list[self.current_scan] = 0
        self.c_list[self.current_scan] = 0
        self.dp50_list[self.current_scan] = (0, 0)
        self.update_scan_info_and_graphs()
        self.move_progress_bar_forward(complete=True)


def draw_all_scans_alignment_summary_graph(self):
    pass
    """
    A graph of peak alignment, and also allow interaction to select peak to process
    """
    # for i in range(len(self.min_pos_CCNC_list)):
    #     if self.min_pos_CCNC_list[i] and self.min_pos_SMPS_list[i] and self.is_usable_for_sigmoid_fit_list[i]:
    #         temp_smps_all_scans_list.append(self.min_pos_SMPS_list[i])
    #         temp_ccnc_all_scans_list.append(self.min_pos_CCNC_list[i])
    #     else:
    #         if len(temp_smps_all_scans_list) > 0:
    #             temp_smps_all_scans_list.append(temp_smps_all_scans_list[-1] + self.scan_duration)
    #             temp_ccnc_all_scans_list.append(temp_ccnc_all_scans_list[-1] + self.scan_duration)
    #         else:
    #             temp_smps_all_scans_list.append(10)
    #             temp_ccnc_all_scans_list.append(10)
    # colors = []
    # for i in range(len(self.is_usable_for_sigmoid_fit_list)):
    #     if self.is_usable_for_sigmoid_fit_list[i]:
    #         colors.append("#43A047")
    #     else:
    #         colors.append("#FF0000")
    # x = numpy.asarray(temp_smps_all_scans_list)
    # y = numpy.asarray(temp_ccnc_all_scans_list)
    #
    # if self.all_scans_alignment_ax is None:
    #     # result = scipy.stats.linregress(x, y)
    #     # slope = result[0]
    #     # yIntercept = result[1]
    #     # Recalculate the position of the smps
    #     # plt.plot(x, x * slope + yIntercept, linewidth=1, color='#43A047', label="Regression line")
    #     # slope = ('{0:.4f}'.format(slope))
    #     # yIntercept = ('{0:.4f}'.format(yIntercept))
    #     # textToShow = str(slope) + "* x" + " + " + str(yIntercept)
    #     # ax.text(10,10,textToShow)
    #     # plt.scatter(x, y, s=150, marker="o", color=colors, picker=6)
    #     figure, ax = plt.subplots(facecolor=settings.GRAPH_BACKGROUND_COLOR)
    #     ax.axes.set_frame_on(False)
    #     ax.grid(color='0.5')
    #     ax.axhline(0, color='0.6', linewidth=4)
    #     ax.axvline(0, color='0.6', linewidth=4)
    #     ax.tick_params(axis='x', color='1', which='both', labelcolor="0.6")
    #     ax.tick_params(axis='y', color='1', which='both', labelcolor="0.6")
    #     ax.yaxis.label.set_color('0.6')
    #     ax.xaxis.label.set_color('0.6')
    #     figure.canvas.mpl_connect('pick_event', self.on_pick)
    #     ax.set_xlabel("SMPS minumum point")
    #     ax.set_ylabel("CCNC minimum point")
    #     handles, labels = ax.get_legend_handles_labels()
    #     legend = ax.legend(handles, labels, loc="upper right", bbox_to_anchor=(1, 0.7))
    #     legend.get_frame().set_facecolor('#9E9E9E')
    #     self.all_scans_alignment_bars = ax.bar(range(1, self.number_of_scan + 1), self.shift_factor_list,
    #                                            color="#43A047", picker=True, align='center')
    #     for i in range(len(self.is_usable_for_sigmoid_fit_list)):
    #         if not self.is_usable_for_sigmoid_fit_list[i]:
    #             self.all_scans_alignment_bars[i].set_facecolor('#111111')
    #     self.all_scans_alignment_figure = figure
    #     self.all_scans_alignment_ax = ax
    #
    # self.view.update_temp_and_min_figure(self.all_scans_alignment_figure)



    # class AverageKappaPointsDataTable(QTableWidget):
# def resize(self, parent_width, parent_height):
    #     self.setFixedHeight(parent_height)
    #     self.setFixedWidth(parent_width)
    #     self.verticalHeader().setDefaultSectionSize(self.height() / 25)

    # def __init__(self, main_window=None):
    #     QTableWidget.__init__(self)
    #     self.setRowCount(0)
    #     self.setShowGrid(False)
    #     self.verticalHeader().setVisible(False)
    #     self.horizontalHeader().setVisible(False)
    #     self.verticalHeader().setDefaultSectionSize(self.height() / 25)
    #     self.setWordWrap(True)
    #     self.setFrameStyle(QFrame.NoFrame)
    #
    #     self.main_window = main_window
    #
    #     # set background color
    #     self.setAutoFillBackground(True)
    #     palette = QPalette()
    #     palette.setColor(QPalette.Base, settings.infoAreaBackgroundColor)
    #     self.setPalette(palette)

    # def update_data(self):
    #     for i in range(self.rowCount()):
    #         self.removeRow(self.rowCount() - 1)
    #
    #     data_dict = self.main_window.controller.alpha_pinene_dict
    #     self.setColumnCount(len(data_dict.keys()) + 1)
    #     self.headerList = []
    #     ssHeader = SingleTableHeaderItem("SS(%)")
    #     meanDPHeader = SingleTableHeaderItem("MeanDp(nm)")
    #     stdDPheader = SingleTableHeaderItem("StdDp(nm)")
    #     meanAppHeader = SingleTableHeaderItem("Mean K, app")
    #     stdAppHeader = SingleTableHeaderItem("Std K, app")
    #     meanAnaHeader = SingleTableHeaderItem("Mean K, ana")
    #     stdAnaHeader = SingleTableHeaderItem("Std K, ana")
    #     meanDeviHeader = SingleTableHeaderItem("Mean of %Deviation")
    #     deviMeanHeader = SingleTableHeaderItem("%Deviation of Mean")
    #     self.headerList.append(ssHeader)
    #     self.headerList.append(meanDPHeader)
    #     self.headerList.append(stdDPheader)
    #     self.headerList.append(meanAppHeader)
    #     self.headerList.append(stdAppHeader)
    #     self.headerList.append(meanAnaHeader)
    #     self.headerList.append(stdAnaHeader)
    #     self.headerList.append(meanDeviHeader)
    #     self.headerList.append(deviMeanHeader)
    #
    #     # Insert Header
    #     for i in range(0, len(self.headerList)):
    #         self.insertRow(self.rowCount())
    #         self.setCellWidget(self.rowCount() - 1, 0, self.headerList[i])
    #     count = 1
    #     for a_key in data_dict.keys():
    #         a_list = [a_key]
    #         a_list.extend(data_dict[a_key])
    #         self.add_message(a_list, count)
    #         count += 1

    # def add_message(self, data_list, column_pos):
    #     for i in range(len(data_list)-1):
    #         a_cell = '% .2f' % data_list[i]
    #         a_cell = SingleTableItem(a_cell)
    #         self.setCellWidget(i, column_pos, a_cell)


class KappaControlTabWidget(QWidget):
    def resize(self, parent_width, parent_height):
        self.setFixedHeight(parent_height * 1 / 10)
        self.setFixedWidth(parent_width)
        self.show_parameters_button.resize(self.width(), self.height())
        self.toggle_average_all_k_points_button.resize(self.width(), self.height())

    def __init__(self, main_window=None):
        self.main_window = main_window
        QWidget.__init__(self)
        self.layout = QHBoxLayout()
        self.layout.setContentsMargins(0, 10, 60, 10)

        self.show_parameters_button = CustomButton("Parameters", main_window)
        self.toggle_average_all_k_points_button = CustomButton("Ave Points", main_window)

        self.show_parameters_button.clicked.connect(self.on_click_show_parameters)
        self.toggle_average_all_k_points_button.clicked.connect(self.on_click_toggle_all_average_k_points)

        self.layout.addWidget(self.show_parameters_button)
        self.layout.addWidget(self.toggle_average_all_k_points_button)

        self.setLayout(self.layout)
        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Background, settings.controlAreaBackgroundColor)
        self.setPalette(palette)
        # self.toggle_all_less_k_lines_button = CustomButton("Less Lines", main_window)
        # self.show_raw_data_button = CustomButton("Kappa Data", main_window)
        # self.show_graph_data_button = CustomButton("Graph Data", main_window)
        # self.show_raw_data_button.clicked.connect(self.on_click_show_raw_data)
        # self.show_graph_data_button.clicked.connect(self.on_click_show_graph)
        # self.layout.addWidget(self.toggle_all_less_k_lines_button)
        # self.layout.addWidget(self.show_raw_data_button)
        # self.layout.addWidget(self.show_graph_data_button)
        # self.toggle_all_less_k_lines_button.clicked.connect(self.on_click_toggle_k_lines)

    def on_click_show_parameters(self):
        self.main_window.centralWidget().info_widget.change_to_parameters_data_table()

    def on_click_toggle_all_average_k_points(self):
        # order of these functions are important. we need to call draw_kappa_graph first
        # to update the list of kappa points, then call show_data_for... to update the UI
        if (self.main_window.controller.is_show_all_k_points):
            self.main_window.controller.is_show_all_k_points = False
            self.toggle_average_all_k_points_button.setText("All Points")
            self.main_window.controller.draw_kappa_graph()
            self.main_window.centralWidget().info_widget.show_data_for_ave_k_points()
        else:
            self.main_window.controller.is_show_all_k_points = True
            self.toggle_average_all_k_points_button.setText("Ave Points")
            self.main_window.controller.draw_kappa_graph()
            self.main_window.centralWidget().info_widget.show_data_for_all_k_points()


    # def on_click_toggle_k_lines(self):
    #     if (self.main_window.controller.is_show_all_k_lines):
    #         self.main_window.controller.is_show_all_k_lines = False
    #         self.toggle_all_less_k_lines_button.setText("All Lines")
    #     else:
    #         self.main_window.controller.is_show_all_k_lines = True
    #         self.toggle_all_less_k_lines_button.setText("Less Lines")
    #     self.main_window.controller.draw_kappa_graph()

    # def on_click_show_graph(self):
    #     self.main_window.centralWidget().info_widget.show_data_for_ave_k_points()

    # def on_click_show_raw_data(self):
    #     self.main_window.centralWidget().info_widget.show_data_for_all_k_points()


# class SingleTableHeaderItem(QWidget):
#     def __init__(self, message):
#         QWidget.__init__(self)
#         self.layout = QHBoxLayout()
#         self.field_text = FieldText(message)
#         self.layout.setContentsMargins(0, 0, 0, 0)
#         self.layout.addWidget(self.field_text)
#         self.setLayout(self.layout)
#         self.setAutoFillBackground(True)
#         palette = QPalette()
#         palette.setColor(QPalette.Base, settings.infoAreaItemBackgroundColor)
#         self.setPalette(palette)
#
#     def toggle_color(self):
#         self.field_text.toggle_color()

class KappaParamsDataTable(QTableWidget):
    def resize(self, parent_width, parent_height):
        self.setFixedHeight(parent_height)
        self.setFixedWidth(parent_width)
        self.verticalHeader().setDefaultSectionSize(self.height() / 25)

    def __init__(self, main_window=None):
        QTableWidget.__init__(self)
        self.setColumnCount(1)
        self.setRowCount(0)
        self.setShowGrid(False)
        self.verticalHeader().setVisible(False)
        self.horizontalHeader().setVisible(False)
        self.horizontalHeader().setStretchLastSection(True)
        self.verticalHeader().setDefaultSectionSize(self.height() / 25)
        self.setWordWrap(True)
        self.resizeRowsToContents()
        self.setFrameStyle(QFrame.NoFrame)
        self.mainWindow = main_window

        self.setAutoFillBackground(True)
        palette = QPalette()
        palette.setColor(QPalette.Base, settings.infoAreaBackgroundColor)
        self.setPalette(palette)

    def update_data(self):
        for i in range(self.rowCount()):
            self.removeRow(self.rowCount() - 1)

        header = AllTableHeader("Variables")
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount() - 1, 0, header)
        controller = self.mainWindow.controller
        kappaVars = (controller.sigma, controller.temp, controller.dd, controller.iKappa, controller.dd2, controller.iKappa2, controller.solubility)
        s = kappaVars[0]
        t = kappaVars[1]
        d1 = kappaVars[2]
        i1 = kappaVars[3]
        d2 = kappaVars[4]
        i2 = kappaVars[5]
        solu = kappaVars[6]
        # TODO: get unit for this
        self.add_message("Sigma", s)
        self.add_message("Temperature (k)", t)
        self.add_message("dry diamater(1) (nm)", d1)
        self.add_message("iKppa(1)", i1)
        self.add_message("dry diamater(2) (nm)", d2)
        self.add_message("iKppa(2)", i2)
        self.add_message("Solubility", solu)

    def add_message(self, field, message):
        message = str(message)
        item = AlignmentTableItem(field, message)
        self.insertRow(self.rowCount())
        self.setCellWidget(self.rowCount()-1, 0, item)

    def toggle_color(self,row):
        pass

    #
    # def switch_to_kappa_widget(self):
    #     self.clear_layout(self.layout)
    #     self.info_widget = KappaInformationAndDataWidget(self.main_window)
    #     self.graph_widget = KappaGraphWidget(self.main_window)
    #     self.layout.addWidget(self.info_widget)
    #     self.layout.addWidget(self.graph_widget)
    #     self.setLayout(self.layout)
    #     self.resize()
    #
    # def switch_to_scan_widget(self):
    #     self.clear_layout(self.layout)
    #     self.info_widget = ScanInformationWidget(self.main_window)
    #     self.graph_widget = ScanGraphsWidget(self.main_window)
    #     self.layout.addWidget(self.info_widget)
    #     self.layout.addWidget(self.graph_widget)
    #     self.setLayout(self.layout)
    #     self.resize()
    #
    # def clear_layout(self, layout):
    #     if layout is not None:
    #         while layout.count():
    #             item = layout.takeAt(0)
    #             widget = item.widget()
    #             if widget is not None:
    #                 widget.deleteLater()
    #             else:
    #                 self.clear_layout(item.layout())
    #
    # def keyReleaseEvent(self, event):
    #     # finish kappa
    #     if self.main_window.controller.kappa_ax:
    #         self.main_window.controller.on_key_release_kappa_graph(event)
    #     else:
    #         self.main_window.controller.on_key_release(event)
#
# class CentralWidgetKappaCalc(QWidget):
#     def __init__(self, parent, raw_conc_time_graph, smoothed_conc_time_graph, ratio_dp_graph, temp_graph):
#         QWidget.__init__(self)
#         # Add widgets
#         # init the necessary contents
#         self.raw_conc_time_graph = raw_conc_time_graph
#         self.smoothed_conc_time_graph = smoothed_conc_time_graph
#         self.ratio_dp_graph = ratio_dp_graph
#         self.temp_graph = temp_graph
#         self.h_splitter_1 = QSplitter(Qt.Horizontal)
#         self.h_splitter_1.addWidget(self.smoothed_conc_time_graph)
#         self.h_splitter_1.addWidget(self.raw_conc_time_graph)
#         self.h_splitter_2 = QSplitter(Qt.Horizontal)
#         self.h_splitter_2.addWidget(self.ratio_dp_graph)
#         self.h_splitter_2.addWidget(self.temp_graph)
#         self.v_splitter = QSplitter(Qt.Vertical)
#         self.v_splitter.addWidget(self.h_splitter_1,self.h_splitter_2)
#
#     def switch_to_kappa_widget(self):
#         self.clear_layout(self.layout)
#         self.info_widget = KappaInformationAndDataWidget(self.main_window)
#         self.graph_widget = KappaGraphWidget(self.main_window)
#         self.layout.addWidget(self.info_widget)
#         self.layout.addWidget(self.graph_widget)
#         self.setLayout(self.layout)
#         self.resize()
#
#     def switch_to_scan_widget(self):
#         self.clear_layout(self.layout)
#         self.info_widget = ScanInformationWidget(self.main_window)
#         self.graph_widget = ScanGraphsWidget(self.main_window)
#         self.layout.addWidget(self.info_widget)
#         self.layout.addWidget(self.graph_widget)
#         self.setLayout(self.layout)
#         self.resize()
#
#     def clear_layout(self, layout):
#         if layout is not None:
#             while layout.count():
#                 item = layout.takeAt(0)
#                 widget = item.widget()
#                 if widget is not None:
#                     widget.deleteLater()
#                 else:
#                     self.clear_layout(item.layout())
#
#     def keyReleaseEvent(self, event):
#         # finish kappa
#         if self.main_window.controller.kappa_ax:
#             self.main_window.controller.on_key_release_kappa_graph(event)
#         else:
#             self.main_window.controller.on_key_release(event)


    # def check_for_update(self):
    #     """
    #     Check for update by accessing a file callede APP_VERSION.html on github. Very much hard-coded
    #     :return:
    #     """
    #     try:
    #         response = urllib2.urlopen(
    #             'https://raw.githubusercontent.com/khaiquangnguyen/Chemics/master/APP_VERSION.html')
    #         html = response.read()
    #         update = int(html[0])
    #         self.show_update_dialog(update)
    #     except urllib2.URLError, e:
    #         # For Python 2.6
    #         if isinstance(e.reason, socket.timeout):
    #             pass
    #         else:
    #             # reraise the original error
    #             raise
    #     except socket.timeout, e:
    #         # For Python 2.7
    #         pass
