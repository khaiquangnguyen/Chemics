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
