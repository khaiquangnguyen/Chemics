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
        self.update_view()
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
        self.update_view()
        self.move_progress_bar_forward(complete=True)


def draw_all_scans_alignment_summary_graph(self):
    pass
    """
    A graph of peak alignment, and also allow interaction to select peak to process
    """
    # for i in range(len(self.min_pos_CCNC_list)):
    #     if self.min_pos_CCNC_list[i] and self.min_pos_SMPS_list[i] and self.usable_for_sigmoid_fit_list[i]:
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
    # for i in range(len(self.usable_for_sigmoid_fit_list)):
    #     if self.usable_for_sigmoid_fit_list[i]:
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
    #     figure, ax = plt.subplots(facecolor=settings.graphBackgroundColor)
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
    #     for i in range(len(self.usable_for_sigmoid_fit_list)):
    #         if not self.usable_for_sigmoid_fit_list[i]:
    #             self.all_scans_alignment_bars[i].set_facecolor('#111111')
    #     self.all_scans_alignment_graph = figure
    #     self.all_scans_alignment_ax = ax
    #
    # self.view.update_temp_and_min_figure(self.all_scans_alignment_graph)
