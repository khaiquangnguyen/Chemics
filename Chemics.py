import re
import pandas
from datetime import *
from scipy import signal
import numpy
from scipy import *
from scipy import stats
from GUI import *
import scipy.constants
import matplotlib.pyplot as plt
import matplotlib
import settings
import scipy.optimize as opt
import time
from PySide import QtGui
from Exceptions import *
import FileDialog
import Tkinter
import copy
import gc
from timeit import default_timer as timer
from HelperFunctions import *
import FastDpCalculator
from settings import *

matplotlib.style.use('ggplot')


class Controller():
    def __init__(self, mode=True):
        """"
        The controller of the program. Also works doubly as the model.
        """

        self.kappa_graph = None
        self.kappa_ax = None
        self.is_show_all_k_lines = False
        self.temperature_graph = None
        self.temperature_axis = None
        self.t1_line = None
        self.t2_line = None
        self.t3_line = None
        self.k_lines_list = []
        self.kappa_points = None
        self.concentration_over_scan_time_graph = None
        self.concentration_over_scan_time_axis = None
        self.scan_time_figure_smps_line = None
        self.scan_time_figure_ccnc_line = None

        self.ccn_cn_ratio_graph = None
        self.ccn_cn_ratio_ax = None
        self.ccn_cn_ratio_corrected_points = None
        self.sigmoid_line = None

        self.sigmoid_fit_graph = None
        self.sigmoid_fit_ax = None
        self.normalized_concentration_points = None
        self.ccn_cn_points = None

        self.min_compare_graph = None

        self.ccnc_sig_list_list = []

        self.view = None
        self.files = None
        self.smps_txt_files = []
        self.ccnc_csv_files = []
        self.current_scan = 0
        self.number_of_scan = 0
        self.scan_start_time_list = None
        self.scan_end_time_list = None
        self.experiment_date = None
        self.scan_duration = 0
        # ccnc = Cloud Condensation Nuclei Counter
        self.ccnc_data = None
        # smps = Scanning Mobility Particle Sizer
        self.smps_data = None
        self.particle_diameter_list = None
        self.normalized_concentration_list = None
        self.super_saturation_list = []
        self.processed_data = None
        self.raw_data = None
        self.current_point = None
        self.finish_scan_alignment = False
        self.cancelling_progress_bar = False
        self.usable_for_sigmoid_fit_list = []

        self.flow_rate = 0.3
        self.ccn_list = None
        self.cn_list = None
        self.ccn_fixed_list = None
        self.cn_fixed_list = None
        self.g_cn_list = None
        self.g_ccn_list = None
        self.ccnc_sig_list = []
        self.ccn_normalized_list = []
        self.diameter_midpoint_list = []
        self.ccn_cn_sim_list = []
        self.temp1 = []
        self.temp2 = []
        self.temp3 = []

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

        self.min_pos_SMPS_list = []
        self.min_pos_CCNC_list = []

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
        self.shift_factor_list = []
        self.invalid_k_points_list = []
        self.kappa_points_data_list = []  # format of data is (dp,ss)
        self.klines_data = None
        self.max_kappa = None
        self.min_kappa = None
        self.is_show_all_k_points = True
        self.kappa_calculate_dict = {}
        self.alpha_pinene_dict = {}
        self.scan_position_in_data = []
        self.usable_for_kappa_cal_list = []
        self.clone = None
        self.all_scans_alignment_graph = None
        self.all_scans_alignment_ax = None
        self.all_scans_alignment_bars = []
        self.all_scans_alignment_visited = []
        self.graph_current_selection = None
        self.kappa_points_is_included_list = {}  # format of key is (dp,ss)

    ##############################################
    #
    # PARSE FILES AND MATCH SMPS AND CCNC DATA
    #
    ##############################################

    def get_concentration(self):
        """
        The flow_rate of particles in the flow.
        In the original Post-Processing Spreadsheet, this is called CountsToConcConv
        :return:
        """
        self.flow_rate = self.view.get_concentration()

    def get_raw_data(self):
        """
        Get the raw SMPS and CCNC data from the files. This converts the data files
        to Python objects and acquire certain necessary information. The following information are acquired
        after this process:
        + self.smps_data: smps data as a python object
        + self.ccnc_data: ccnc data as a python object
        + self.experiment_date: the date of the experiment from which we have the data
        + self.scan_duration: the duration of each scan
        + self.number_of_scan: the total number of scans
        + self.scan_start_time_list: the list of the start time of each scan
        + self.scan_end_time_list: the list of the end time of each scan
        """
        self.smps_txt_files = []
        self.ccnc_csv_files = []
        for a_file in self.files:
            if a_file.lower().endswith('.txt'):
                self.smps_txt_files.append(a_file)
            elif a_file.lower().endswith('.csv'):
                self.ccnc_csv_files.append(a_file)
        if (len(self.smps_txt_files) == 0) or (len(self.ccnc_csv_files) == 0):
            raise FileNotFoundError()
        else:
            self.smps_txt_files = [str(x) for x in self.smps_txt_files]
            self.smps_txt_files = self.smps_txt_files[0]
            self.ccnc_csv_files = [str(x) for x in self.ccnc_csv_files]

        try:
            self.ccnc_data = process_csv_files(self.ccnc_csv_files)
            self.smps_data = process_text_files(self.smps_txt_files)
            self.experiment_date = self.ccnc_data[0]
            self.scan_start_time_list = self.smps_data[0]
            for i in range(len(self.smps_data)):
                if ''.join(self.smps_data[i][0].split()).lower() == "scanuptime(s)":
                    scan_up_time = self.smps_data[i][1]
                    scan_down_time = self.smps_data[i + 1][1]
                    break
            self.scan_duration = int(scan_up_time) + int(scan_down_time)
            self.scan_end_time_list = []
            for i in range(len(self.scan_start_time_list)):
                self.scan_end_time_list.append(datetime.strptime(self.scan_start_time_list[i], "%H:%M:%S") + timedelta(
                    seconds=self.scan_duration))
            # Remove the scans with invalid start time and end time. This is often due to the mismatch between
            # SMPS and CCNC data.
            while (True):
                smps_end_time = self.scan_end_time_list[-1]
                ccnc_end_time = datetime.strptime(self.ccnc_data[1][len(self.ccnc_data[1]) - 1][0], "%H:%M:%S")
                if smps_end_time > ccnc_end_time:
                    self.scan_start_time_list = self.scan_start_time_list[:-1]
                    self.scan_end_time_list = self.scan_end_time_list[:-1]
                else:
                    break
            self.scan_end_time_list = [datetime.strftime(x, "%H:%M:%S") for x in self.scan_end_time_list]
            self.number_of_scan = len(self.scan_start_time_list)
        except:
            raise FileParseError()

    def get_normalized_concentration(self):
        """
        Get the normalized flow_rate data. This information is taken directly from the SMPS (.txt) file
        Normalized flow_rate is also called dN/dLogDp, which is the notation we use for the graph.
        For better understanding, you can read the following document.
        http://www.tsi.com/uploadedFiles/Product_Information/Literature/Application_Notes/PR-001-RevA_Aerosol-Statistics-AppNote.pdf
        :return:
        """
        try:
            for k in range(1, len(self.smps_data)):
                if re.search('[aParam-zA-Z]', self.smps_data[k][0]):
                    pos_in_file = k
                    break
            startTime = ['dp'] + self.scan_start_time_list
            endTime = ['dp'] + self.scan_end_time_list
            self.normalized_concentration_list = [startTime, endTime]
            for k in range(1, pos_in_file):
                self.normalized_concentration_list.append(self.smps_data[k][:len(startTime)])
        except:
            raise NormalizedConcentrationError()

    def get_smps_count_and_ccnc_count(self):
        """
        Get the SMPS count from the SMPS files and CCNC count from CCNC files.
        The SMPS count is taken directly from the SMPS files, but the CCNC count requires some processing
        on the data from the CCNC file.
        :return:
        """
        ccnc_data = self.ccnc_data[1]
        smps_data = self.smps_data
        startTime = self.scan_start_time_list
        endTime = self.scan_end_time_list
        width = self.number_of_scan
        try:
            # get the position (line number) of the SMPS count information in the SMPS file
            for k in range(3, len(self.smps_data)):
                if re.search('[aParam-zA-Z]', self.smps_data[k][0]):
                    for j in range(k, len(self.smps_data)):
                        if not re.search('[aParam-zA-Z]', self.smps_data[j][0]):
                            countPos = j
                            break
                    break
            # Get the SMPS count
            loopCount = 0
            countList = []
            smpsList = []
            sumSMPS = 0
            count = 0
            for k in range(countPos, len(smps_data) - 1):
                loopCount += 1
                for j in range(0, width):
                    sumSMPS += float(smps_data[k][j * 2 + 1]) * int(smps_data[k][j * 2 + 2])
                if not countList:
                    for j in range(1, width + 1):
                        num = int(smps_data[k][j * 2])
                        countList.append(num)
                else:
                    for j in range(1, width + 1):
                        num = int(smps_data[k][j * 2])
                        countList[j - 1] += num
                if loopCount == 10:
                    count += 1
                    if sum(countList) == 0:
                        smpsList.append([float(smps_data[k][1])] + [count] + countList)
                    else:
                        smpsList.append([sumSMPS / sum(countList)] + [count] + countList)
                    loopCount = 0
                    sumSMPS = 0
                    countList = []
        except:
            raise SMPSCountError()

        try:
            # Steps the get the CCNC count
            ccnList = []
            aveSizeList = []
            extraCCNList = []
            extraAveSizeList = []
            sizeList = [0.625] + [0.875]
            size = 1.25
            binPos = 25
            aSSList = []
            # Calculate the size of each bin
            for i in range(0, 18):
                sizeList.append(size)
                size += 0.5
            # Get the first position of CCNC count in the ccnc file
            while (True):
                timeStamp = startTime[0]
                if ccnc_data[0][0] > timeStamp:
                    startTime = startTime[1:]
                    endTime = endTime[1:]
                else:
                    break
            for k in range(0, len(ccnc_data) - 1):
                if ccnc_data[k][0] == timeStamp and ccnc_data[k + 1][0] != timeStamp:
                    break

            i = 0
            t = 0
            isPeakTimeStamp = True
            smpsCcnList = []
            timeStamp = datetime.strptime(startTime[0], "%H:%M:%S")
            self.scan_position_in_data = []
            # Get the CCNC count
            while (True):
                sizeSum = 0
                countSum = 0
                previousTimeStamp = datetime.strptime(startTime[i], "%H:%M:%S").time()
                super_saturation = -1
                # Get all data within the duration of a scan
                self.scan_position_in_data.append(len(smpsCcnList))
                for t in range(self.scan_duration):
                    aLine = [smpsList[t][1]] + [timeStamp.time()] + [float(smpsList[t][0])] + [
                        float(smpsList[t][i + 2])]
                    sizeSum = 0
                    countSum = 0
                    aLine.append(float(ccnc_data[k + t][-3]))
                    if t == 0:
                        super_saturation = float(ccnc_data[k + t][1])
                    elif float(ccnc_data[k + t][1]) != super_saturation:
                        super_saturation = -1
                    for m in range(0, 20):
                        sizeSum += sizeList[m] * float(ccnc_data[k + t][binPos + m])
                        size += 0.5
                        countSum += float(ccnc_data[k + t][binPos + m])
                        # get the average size for each scan
                    if countSum == 0:
                        aLine.append(0)
                    else:
                        aLine.append(sizeSum / countSum)
                    aLine.append(float(ccnc_data[k + t][5]))
                    aLine.append(float(ccnc_data[k + t][7]))
                    aLine.append(float(ccnc_data[k + t][9]))
                    smpsCcnList.append(aLine)
                    timeStamp += timedelta(seconds=1)
                k += self.scan_duration

                # Update the lists
                self.super_saturation_list.append(super_saturation)

                # If reach the end of smps, collect the rest of ccnc count
                # so that we can use them latter for scan alignment
                timeGap = 0
                if i == len(startTime) - 1:
                    timeGap = (datetime.strptime(ccnc_data[-1][0], "%H:%M:%S") - timeStamp).seconds
                # Do whatever
                else:
                    nextTimeStamp = datetime.strptime(startTime[i + 1], "%H:%M:%S")
                    if nextTimeStamp < timeStamp:
                        timeGap = -(timeStamp - nextTimeStamp).seconds
                    else:
                        timeGap = (nextTimeStamp - timeStamp).seconds

                # If the time shown in header is smaller than the supposed actual time frame
                # Then correct the header accordingly
                if timeGap < 0:
                    startTime[i + 1] = (
                        datetime.strptime(startTime[i + 1], "%H:%M:%S") + timedelta(seconds=abs(timeGap))).strftime(
                        "%I:%M:%S")
                    self.scan_start_time_list = startTime
                    timeGamp = 0
                for t in range(timeGap):
                    aLine = [None] + [timeStamp.time()] + [None] + [None]
                    sizeSum = 0
                    countSum = 0
                    aLine.append(float(ccnc_data[k + t][-3]))
                    if ccnc_data[k + t][1] > super_saturation:
                        super_saturation = ccnc_data[k + t][1]
                    for m in range(0, 20):
                        sizeSum += sizeList[m] * float(ccnc_data[k + t][binPos + m])
                        size += 0.5
                        countSum += float(ccnc_data[k + t][binPos + m])
                        # get the average size for each scan
                    if countSum == 0:
                        aLine.append(0)
                    else:
                        aLine.append(sizeSum / countSum)
                    aLine.append(float(ccnc_data[k + t][5]))
                    aLine.append(float(ccnc_data[k + t][7]))
                    aLine.append(float(ccnc_data[k + t][9]))
                    smpsCcnList.append(aLine)
                    timeStamp += timedelta(seconds=1)
                k += timeGap
                # Loop break condition
                i += 1
                if i >= len(startTime):
                    break
        except:
            raise CCNCCountError()

        title = ["scan time"] + ["real time"] + ["dp"] + ["SMPS Count"] + \
                ["CCNC Count"] + ["Ave Size"] + ["T1"] + ["T2"] + ["T3"]
        self.raw_data = pandas.DataFrame(smpsCcnList, columns=title)

    def match_smps_ccnc_data(self):
        """
        The most important step of the whole program. This step will detects the mismatch of time between the SMPS
        and the CCNC data, and will correct this by moving the CCNC count backward or forward a certain number of
        seconds. The SMPS and CCNC data are considered matched when their local minimum between the two peaks of a
        scan are at the same time location.
        :return:
        """
        self.usable_for_kappa_cal_list = []
        self.usable_for_sigmoid_fit_list = []
        self.min_pos_CCNC_list = []
        self.min_pos_SMPS_list = []
        num_additional_data_point = int(self.scan_duration / 2)
        start_time = 0
        end_time = self.scan_duration
        current_scan = 0
        new_data = [self.raw_data.columns.values.tolist()]
        min_dist = self.scan_duration / 10
        shift_factor = 0
        self.shift_factor_list = []
        self.number_of_scan = 0
        while True:
            a_scan = numpy.asarray(self.raw_data[start_time:end_time]["SMPS Count"])
            min_pos_smps = get_min_index(a_scan)
            if min_pos_smps == -1:
                self.min_pos_SMPS_list.append(None)
                self.min_pos_CCNC_list.append(None)
                self.usable_for_sigmoid_fit_list.append(False)
                self.usable_for_kappa_cal_list.append(False)
                min_pos_smps = 0
                min_pos_ccnc = 0
            else:
                self.min_pos_SMPS_list.append(min_pos_smps + start_time)
                a_scan = numpy.asarray(
                    self.raw_data[start_time + shift_factor:end_time + num_additional_data_point]["CCNC Count"])
                min_pos_ccnc = get_min_index(a_scan)
                if min_pos_ccnc == -1:
                    self.min_pos_SMPS_list[-1] = None
                    self.min_pos_CCNC_list.append(None)
                    self.usable_for_sigmoid_fit_list.append(False)
                    self.usable_for_kappa_cal_list.append(False)
                    min_pos_smps = 0
                    min_pos_ccnc = 0
                else:
                    self.usable_for_kappa_cal_list.append(True)
                    self.usable_for_sigmoid_fit_list.append(True)
                    shift_factor = shift_factor + min_pos_ccnc - min_pos_smps
                    if current_scan == 0:
                        self.min_pos_CCNC_list.append(min_pos_smps)
                    else:
                        self.min_pos_CCNC_list.append(min_pos_ccnc + start_time + (min_pos_ccnc - min_pos_smps))
            self.shift_factor_list.append(shift_factor)

            for i in range(self.scan_duration):
                scan_time = i + 1
                real_time = self.raw_data.iat[start_time + i, 1]
                dp = self.raw_data.iat[start_time + i, 2]
                smpsCount = self.raw_data.iat[start_time + i, 3]
                ccncCount = self.raw_data.iat[start_time + i + shift_factor, 4]
                aveSize = self.raw_data.iat[start_time + i + shift_factor, 5]
                t1 = self.raw_data.iat[start_time + i + shift_factor, 6]
                t2 = self.raw_data.iat[start_time + i + shift_factor, 7]
                t3 = self.raw_data.iat[start_time + i + shift_factor, 8]
                new_data.append([scan_time, real_time, dp, smpsCount, ccncCount, aveSize, t1, t2, t3])
            current_scan += 1
            self.number_of_scan += 1
            if end_time + self.scan_duration + shift_factor >= len(self.raw_data) or current_scan >= len(
                    self.scan_position_in_data):
                break

            start_time = self.scan_position_in_data[current_scan]
            end_time = start_time + self.scan_duration
        headers = new_data.pop(0)
        self.processed_data = pandas.DataFrame(new_data, columns=headers)
        self.processed_data = convert_date(self.processed_data)

    def parse_and_match_smps_ccnc_data(self):
        try:
            self.reset()
            self.view.reset()
            self.move_progress_bar_forward(max_value=5)
            self.get_concentration()
            self.move_progress_bar_forward("Processing SMPS and CCNC files...")
            self.get_raw_data()
            self.move_progress_bar_forward("Acquiring normalized flow_rate...")
            self.get_normalized_concentration()
            self.move_progress_bar_forward("Acquiring SMPS count and CCNC count...")
            self.get_smps_count_and_ccnc_count()
            self.move_progress_bar_forward("Matching SMPS and CCNC data....")
            self.match_smps_ccnc_data()
            self.move_progress_bar_forward(complete=True)
            self.move_progress_bar_forward(max_value=self.number_of_scan)
            for i in range(0, self.number_of_scan):
                self.move_progress_bar_forward("Processing scan number " + str(i + 1), value=0)
                plt.ioff()
                self.current_scan = i
                self.prepare_scan_data()
                if (self.super_saturation_list[i] == -1 or not check_temperature_fluctuation(self.temp1)
                    or not check_temperature_fluctuation(self.temp2) or not check_temperature_fluctuation(self.temp3)):
                    self.min_pos_CCNC_list[i] = None
                    self.min_pos_SMPS_list[i] = None
                    self.usable_for_sigmoid_fit_list[i] = False
            self.move_progress_bar_forward(complete=True)
            self.current_scan = -1
            self.view.update_experiment_information()
            self.switch_to_scan(0)
        except FileNotFoundError:
            self.view.show_error_dialog("Fail to find SMPS or CCNC data in files: " + str(self.files))
        except FileParseError:
            self.view.show_error_dialog("Fail to process the SMPS or CCNC files! Please check the files again!")
        except NormalizedConcentrationError:
            self.view.show_error_dialog("Fail to acquire normalized concentration (dN/dLogDp) "
                                        "data from the SMPS(.txt) file!")
        except SMPSCountError():
            self.view.show_error_dialog(
                "Fail to acquire SMPS count from SMPS file! Please check the SMPS(.txt) file again!")
        except CCNCCountError():
            self.view.show_error_dialog(
                "Fail to acquire CCNC count from CCNC file! Please check the CCNC(.csv) files again!")
        except ProgressBarInterruptException:
            self.reset()
            self.view.reset()
            self.view.show_error_dialog("Cancelling parsing SMPS and CCNC files!")

    # ----------------------graphs------------------------
    def draw_temperature_graph(self):
        if self.temperature_axis is None:
            figure, ax = plt.subplots(facecolor=settings.GRAPH_BACKGROUND_COLOR)
            ax.axes.set_frame_on(False)
            ax.grid(color=GRID_COLOR)
            ax.axhline(0, color=AX_LINE_COLOR, linewidth=4)
            ax.axvline(0, color=AX_LINE_COLOR, linewidth=4)
            ax.tick_params(color=AX_TICK_COLOR, which='both', labelcolor=AX_TICK_COLOR, labelsize=AX_TICK_SIZE)
            ax.set_xlabel("Scan time(s)", color=LABEL_COLOR, size=LABEL_SIZE)
            ax.set_ylabel("Temperature", color=LABEL_COLOR, size=LABEL_SIZE)
            ax.set_title("Temperature over scan time", color=TITLE_COLOR, size=TITLE_SIZE)
            x = range(self.scan_duration)
            minY = min(min(self.temp1), min(self.temp2), min(self.temp3)) - 3
            maxY = max(max(self.temp1), max(self.temp2), max(self.temp3)) + 3
            ax.axes.set_ylim([minY, maxY])
            self.t1_line, = ax.plot(x, self.temp1, linewidth=5, color=T1_LINE_COLOR, label="T1")
            self.t2_line, = ax.plot(x, self.temp2, linewidth=5, color=T2_LINE_COLOR, label="T2")
            self.t3_line, = ax.plot(x, self.temp3, linewidth=5, color=T3_LINE_COLOR, label="T3")
            handles, labels = ax.get_legend_handles_labels()
            legend = ax.legend(handles, labels, facecolor=LEGEND_BG_COLOR, fontsize=LEGEND_FONT_SIZE)
            self.temperature_graph = figure
            self.temperature_axis = ax
        else:
            minY = min(min(self.temp1), min(self.temp2), min(self.temp3)) - 3
            maxY = max(max(self.temp1), max(self.temp2), max(self.temp3)) + 3
            self.temperature_axis.axes.set_ylim([minY, maxY])
            self.t1_line.set_ydata(self.temp1)
            self.t2_line.set_ydata(self.temp2)
            self.t3_line.set_ydata(self.temp3)
        self.view.update_temp_and_min_figure(self.temperature_graph)

    def draw_concentration_over_scan_time_graph(self):
        if len(self.cn_list) != len(self.ccn_list) or len(self.cn_list) == 0 or len(self.ccn_list) == 0:
            return
        if self.concentration_over_scan_time_axis is None:
            figure, ax = plt.subplots(facecolor=settings.GRAPH_BACKGROUND_COLOR)
            ax.axes.set_frame_on(False)
            ax.grid(color=GRID_COLOR)
            ax.axhline(0, color=AX_LINE_COLOR, linewidth=4)
            ax.axvline(0, color=AX_LINE_COLOR, linewidth=4)
            ax.tick_params(color=AX_TICK_COLOR, which='both', labelcolor=AX_TICK_COLOR, labelsize=AX_TICK_SIZE)
            ax.set_xlabel("Scan time(s)", color=LABEL_COLOR, size=LABEL_SIZE)
            ax.set_ylabel("Concentration (1/cm3)", color=LABEL_COLOR, size=LABEL_SIZE)
            x = range(self.scan_duration)
            self.scan_time_figure_smps_line, = ax.plot(x, self.cn_list, linewidth=4, color=SMPS_LINE_COLOR,
                                                       label="SMPS")
            self.scan_time_figure_ccnc_line, = ax.plot(x, self.ccn_list, linewidth=4, color=CCNC_LINE_COLOR,
                                                       label="CCNC")
            handles, labels = ax.get_legend_handles_labels()
            legend = ax.legend(handles, labels, facecolor=LEGEND_BG_COLOR, fontsize=LEGEND_FONT_SIZE)
            ax.set_title("SMPS and CCNC concentration over scan time", color=TITLE_COLOR, size=TITLE_SIZE)
            self.concentration_over_scan_time_axis = ax
            self.concentration_over_scan_time_graph = figure
        else:
            self.scan_time_figure_smps_line.set_ydata(self.cn_list)
            self.scan_time_figure_ccnc_line.set_ydata(self.ccn_list)
        self.view.update_alignment_and_sigmoid_fit_figures(self.concentration_over_scan_time_graph,
                                                           None)

    def draw_ccn_cn_ratio_over_diameter_graph(self):
        if self.ccn_cn_ratio_ax is None:
            figure, ax = plt.subplots(facecolor=settings.GRAPH_BACKGROUND_COLOR)
            ax.axes.set_frame_on(False)
            ax.grid(color='0.5')
            ax.axhline(0, color=AX_LINE_COLOR, linewidth=2)
            ax.axvline(0, color=AX_LINE_COLOR, linewidth=4)
            ax.axhline(1, color=str(float(AX_LINE_COLOR) + 0.1), linewidth=2, linestyle='dashed')
            ax.tick_params(color=AX_TICK_COLOR, which='both', labelcolor=AX_TICK_COLOR, labelsize=AX_TICK_SIZE)
            ax.set_xlabel("Diameter (nm)", color=LABEL_COLOR, size=LABEL_SIZE)
            ax.set_ylabel("CCNC/SMPS", color=LABEL_COLOR, size=LABEL_SIZE)
            ax.set_title("CCNC/SMPS over Dry Diameter", color=TITLE_COLOR, size=TITLE_SIZE)
            self.ccn_cn_ratio_points, = ax.plot(self.particle_diameter_list, self.ccn_cn_ratio_list, 'o',
                                                color=CCNC_SMPS_POINT_COLOR, mew=0.5,
                                                ms=9, label="CCNC/SMPS")
            yLim = min(2, max(self.ccn_cn_ratio_list)) + 0.2
            ax.axes.set_ylim([-0.1, yLim])
            handles, labels = ax.get_legend_handles_labels()
            legend = ax.legend(handles, labels, facecolor=LEGEND_BG_COLOR, fontsize=LEGEND_FONT_SIZE)
            self.ccn_cn_ratio_ax = ax
            self.ccn_cn_ratio_graph = figure
        else:
            yLim = min(2, max(self.ccn_cn_ratio_list)) + 0.2
            self.ccn_cn_ratio_ax.axes.set_ylim([-0.1, yLim])
            self.ccn_cn_ratio_points.set_xdata(self.particle_diameter_list)
            self.ccn_cn_ratio_points.set_ydata(self.ccn_cn_ratio_list)
        self.view.update_alignment_and_sigmoid_fit_figures(None,
                                                           self.ccn_cn_ratio_graph)

    ##############################################
    #
    # CORRECT CHARGES AND FIT SIGMOID
    #
    ##############################################

    def init_correct_charges(self):
        """
        Initiate the correct charge procedure
        """
        self.cn_fixed_list = self.cn_list[:]
        self.ccn_fixed_list = self.ccn_list[:]
        self.g_ccn_list = self.ccn_fixed_list[:]
        self.g_cn_list = self.cn_fixed_list[:]

    def get_parameters_for_sigmoid_fit(self, min_dp=None, min_dp_asym=None, max_dp_asym=None):
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
            self.usable_for_kappa_cal_list[self.current_scan] = True
        except:
            raise SigmoidFitGetParameterError()

    def fit_sigmoid_line(self):
        # TODO: show why some runs fail
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
            self.usable_for_kappa_cal_list[self.current_scan] = True
        except:
            raise SigmoidLineFitError()

    def refitting_sigmoid_line(self, min_dry_diameter, min_dry_diameter_asymptote, max_dry_diameter_asymptote):
        self.move_progress_bar_forward("Refitting sigmoid line to scan #" + str(self.current_scan + 1),
                                       max_value=2)
        try:
            if not self.usable_for_sigmoid_fit_list[self.current_scan]:
                self.view.show_error_dialog("Can't fit sigmoid line to current scan!")
                raise SigmoidLineFitError()
            self.prepare_scan_data()
            self.move_progress_bar_forward()
            self.get_parameters_for_sigmoid_fit(min_dry_diameter, min_dry_diameter_asymptote,
                                                max_dry_diameter_asymptote)
            self.fit_sigmoid_line()
            self.move_progress_bar_forward()
        except:
            self.usable_for_kappa_cal_list[self.current_scan] = False
            self.usable_for_sigmoid_fit_list[self.current_scan] = False
            self.b_list[self.current_scan] = 0
            self.d_list[self.current_scan] = 0
            self.c_list[self.current_scan] = 0
            self.dp50_list[self.current_scan] = (0, 0)
            self.dp50_wet_list[self.current_scan] = 0
            self.dp50_plus_20_list[self.current_scan] = 0
            self.min_dp_list[self.current_scan] = 0
            self.min_dp_asym_list[self.current_scan] = 0
            self.max_dp_asym_list[self.current_scan] = 0
            self.update_view()
            self.move_progress_bar_forward(complete=True)
        else:
            self.usable_for_kappa_cal_list[self.current_scan] = True
            self.min_dp_list[self.current_scan] = self.min_dp
            self.min_dp_asym_list[self.current_scan] = self.min_dp_asym
            self.max_dp_asym_list[self.current_scan] = self.max_dp_asym
            self.b_list[self.current_scan] = self.b
            self.d_list[self.current_scan] = self.d
            self.c_list[self.current_scan] = self.c
            self.dp50_list[self.current_scan] = (self.d, self.super_saturation_list[self.current_scan])
            self.dp50_wet_list[self.current_scan] = self.dp50_wet
            self.dp50_plus_20_list[self.current_scan] = self.dp50_plus_20
            self.update_view()
            self.move_progress_bar_forward(complete=True)

    def correct_charges_and_fit_sigmoid_one_scan(self):
        """
        optimize a single run
        """
        try:
            if not self.usable_for_sigmoid_fit_list[self.current_scan]:
                raise SigmoidLineFitError
            self.ccnc_sig_list = []
            self.prepare_scan_data()
            remove_small_ccn(self.ccn_list, self.min_ccn)
            self.particle_diameter_list = numpy.asarray(self.particle_diameter_list)
            self.cn_list = numpy.asarray(self.cn_list)
            self.ccn_list = numpy.asarray(self.ccn_list)
            self.cn_fixed_list = numpy.asarray(self.cn_fixed_list)
            self.ccn_fixed_list = numpy.asarray(self.ccn_fixed_list)
            self.g_cn_list = numpy.asarray(self.g_cn_list)
            self.g_ccn_list = numpy.asarray(self.g_ccn_list)
            self.init_correct_charges()
            for i in range(5):
                self.particle_diameter_list, self.cn_list, self.ccn_list, self.cn_fixed_list, \
                self.ccn_fixed_list, self.g_cn_list, self.g_ccn_list = FastDpCalculator.correct_charges(
                    self.particle_diameter_list, self.cn_list, self.ccn_list, self.cn_fixed_list,
                    self.ccn_fixed_list,
                    self.g_cn_list, self.g_ccn_list)
            for i in range(len(self.ccn_fixed_list)):
                self.ccnc_sig_list.append(self.ccn_fixed_list[i] / self.cn_fixed_list[i])
            self.get_parameters_for_sigmoid_fit()
            self.fit_sigmoid_line()
        except:
            self.usable_for_kappa_cal_list[self.current_scan] = False
            self.b_list.append(0)
            self.d_list.append(0)
            self.c_list.append(0)
            self.dp50_list.append((0, 0))
            self.dp50_wet_list.append(0)
            self.dp50_plus_20_list.append(0)
            self.min_dp_list.append(0)
            self.min_dp_asym_list.append(0)
            self.max_dp_asym_list.append(0)
            self.ccnc_sig_list_list.append(None)
        else:
            # Store processed_data
            self.usable_for_kappa_cal_list[self.current_scan] = True
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
            self.finish_scan_alignment = True
            self.switch_to_scan(0)
        except ProgressBarInterruptException:
            self.view.show_error_dialog("The sigmoid fitting process is cancelled!")
            self.finish_scan_alignment = False
            self.parse_and_match_smps_ccnc_data()

    # -------------------------graphs-----------------------------

    def draw_complete_sigmoid_graph(self, new_figure=None):
        """
        Make complete graph of the dry diameter after optimization and sigmodal fit
        """
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
            yLim = min(2, max(self.ccnc_sig_list)) + 0.2
            ax.axes.set_ylim([-0.1, yLim])
            self.normalized_concentration_points, = plt.plot(self.diameter_midpoint_list, self.ccn_normalized_list,
                                                             linewidth=4, color=NORMALIZED_CONCENTRATION_POINT_COLOR,
                                                             label="dN/dLogDp")
            self.ccn_cn_ratio_points, = ax.plot(self.particle_diameter_list, self.ccn_cn_ratio_list, 'o',
                                                color=CCNC_SMPS_POINT_COLOR,
                                                mew=0.5, ms=9, label="CCNC/SMPS")
            self.ccn_cn_ratio_corrected_points, = ax.plot(self.particle_diameter_list, self.ccnc_sig_list, 'o',
                                                          color=CCNC_SMPS_RATIO_CORRECTED_POINT_COLOR, mew=0.5,
                                                          mec="#0D47A1",
                                                          ms=9, label="CCNC/SMPS (Corrected)")
            if self.usable_for_kappa_cal_list[self.current_scan] and \
                    self.usable_for_sigmoid_fit_list[self.current_scan]:
                self.sigmoid_line, = ax.plot(self.particle_diameter_list, self.ccn_cn_sim_list, linewidth=5,
                                             color=SIGMOID_LINE_COLOR, label="Sigmodal Fit")
            handles, labels = ax.get_legend_handles_labels()
            legend = ax.legend(handles, labels, fontsize=LEGEND_FONT_SIZE, facecolor=LEGEND_BG_COLOR)
            self.sigmoid_fit_graph = figure
            self.sigmoid_fit_ax = ax
        else:
            yLim = min(2, max(self.ccnc_sig_list)) + 0.2
            self.sigmoid_fit_ax.axes.set_ylim([-0.1, yLim])
            self.normalized_concentration_points.set_xdata(self.diameter_midpoint_list)
            self.normalized_concentration_points.set_ydata(self.ccn_normalized_list)
            self.ccn_cn_ratio_points.set_xdata(self.particle_diameter_list)
            self.ccn_cn_ratio_points.set_ydata(self.ccn_cn_ratio_list)
            self.ccn_cn_ratio_corrected_points.set_xdata(self.particle_diameter_list)
            self.ccn_cn_ratio_corrected_points.set_ydata(self.ccnc_sig_list)
            if self.usable_for_kappa_cal_list[self.current_scan] and \
                    self.usable_for_sigmoid_fit_list[self.current_scan]:
                if self.sigmoid_line is not None:
                    self.sigmoid_line.set_xdata(self.particle_diameter_list)
                    self.sigmoid_line.set_ydata(self.ccn_cn_sim_list)
                else:
                    self.sigmoid_line, = self.sigmoid_fit_ax.plot(self.particle_diameter_list, self.ccn_cn_sim_list,
                                                                  linewidth=5, color='#EF5350', label="Sigmodal Fit")
            else:
                if self.sigmoid_line:
                    self.sigmoid_line.set_xdata([])
                    self.sigmoid_line.set_ydata([])

        self.view.update_alignment_and_sigmoid_fit_figures(None, self.sigmoid_fit_graph)

    def draw_all_scans_alignment_summary_graph(self):
        if self.all_scans_alignment_ax is None:
            for i in range(self.number_of_scan):
                self.all_scans_alignment_visited.append(False)

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
            handles, labels = ax.get_legend_handles_labels()
            legend = ax.legend(handles, labels, fontsize=LEGEND_FONT_SIZE, facecolor=LEGEND_BG_COLOR)
            self.all_scans_alignment_bars = ax.bar(range(1, self.number_of_scan + 1), self.shift_factor_list,
                                                   color=SCAN_SUMMARY_USABLE_SCAN_COLOR, picker=True, align='center')
            for i in range(len(self.usable_for_sigmoid_fit_list)):
                if not self.usable_for_sigmoid_fit_list[i] or not self.usable_for_kappa_cal_list[i]:
                    self.all_scans_alignment_bars[i].set_facecolor(SCAN_SUMMARY_UNUSABLE_SCAN_COLOR)
            self.all_scans_alignment_visited[self.current_scan] = True
            for i in range(len(self.all_scans_alignment_visited)):
                if self.all_scans_alignment_visited[i]:
                    self.all_scans_alignment_bars[i].set_alpha(1)
                else:
                    self.all_scans_alignment_bars[i].set_alpha(0.3)
            self.all_scans_alignment_bars[self.current_scan].set_facecolor(SCAN_SUMMARY_HIGHLIGHT_COLOR)
            self.all_scans_alignment_graph = figure
            self.all_scans_alignment_ax = ax
        else:
            for i in range(len(self.usable_for_sigmoid_fit_list)):
                if not self.usable_for_sigmoid_fit_list[i] or not self.usable_for_kappa_cal_list[i]:
                    self.all_scans_alignment_bars[i].set_facecolor(SCAN_SUMMARY_UNUSABLE_SCAN_COLOR)
                else:
                    self.all_scans_alignment_bars[i].set_facecolor(SCAN_SUMMARY_USABLE_SCAN_COLOR)
            self.all_scans_alignment_visited[self.current_scan] = True
            for i in range(len(self.all_scans_alignment_visited)):
                if self.all_scans_alignment_visited[i]:
                    self.all_scans_alignment_bars[i].set_alpha(1)
                else:
                    self.all_scans_alignment_bars[i].set_alpha(0.3)
            self.all_scans_alignment_bars[self.current_scan].set_facecolor(SCAN_SUMMARY_HIGHLIGHT_COLOR)

        self.view.update_temp_and_min_figure(self.all_scans_alignment_graph)

    def on_pick(self, event):
        for i in range(len(self.all_scans_alignment_bars)):
            if event.artist == self.all_scans_alignment_bars[i]:
                if i != self.current_scan:
                    self.switch_to_scan(i)

    def on_key_release(self, event):
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
        if dd:
            self.dd = dd
        if i_kappa:
            self.iKappa = i_kappa
        if dd2:
            self.dd2 = dd2
        if iKapp2:
            self.iKappa2 = i_kappa_2
        if solu:
            self.solubility = solu

    def calculate_all_kappa_values(self):
        """
        Calculate the kappa values for every scan.
        """
        # TODO: export the kappa data to excel
        if self.kappa_excel is None:
            self.kappa_excel = pandas.read_csv("kCal.csv", header=None)
        lookup = self.kappa_excel
        self.aParam = 0.00000869251 * self.sigma / self.temp
        self.asc = (exp(sqrt(4 * self.aParam ** 3 / (27 * self.iKappa * (self.dd * 0.000000001) ** 3))) - 1) * 100
        # Calculate each kappa
        firstAKappa = 0
        scCalcs = False
        for i in range(len(self.dp50_list)):
            if not self.usable_for_kappa_cal_list[i]:
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

    def draw_kappa_graph(self):
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
        if self.current_point is None:
            self.current_point = len(self.kappa_points_data_list) - 1
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
            x = k_points_x_list[self.current_point]
            y = k_points_y_list[self.current_point]
            # graph the current selection point
            self.graph_current_selection, = ax.plot(x, y, 'o', color=KAPPA_CURRENT_SELECTION_COLOR,
                                                    mew=0.5, ms=18, label="current selection")
            handles, labels = ax.get_legend_handles_labels()
            legend = ax.legend(handles, labels, facecolor=LEGEND_BG_COLOR, fontsize=LEGEND_FONT_SIZE)
            self.kappa_graph = figure
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
            x = k_points_x_list[self.current_point]
            y = k_points_y_list[self.current_point]
            self.graph_current_selection.set_xdata(x)
            self.graph_current_selection.set_ydata(y)

    def on_key_release_kappa_graph(self, event):
        """
        When keys are pressed. Used to select current kappa points
        :param event: key press event
        :return:
        """
        key = event.key()
        if self.current_point is None:
            self.current_point = len(self.kappa_points_data_list) - 1
        # left arrow key
        elif (key == 16777234):
            self.current_point = min(len(self.kappa_points_data_list) - 1, self.current_point + 1)
        # right arrow key
        elif (key == 16777236):
            self.current_point = max(0, self.current_point - 1)
        # up arrow key
        elif (key == 16777235):
            curr_ss = self.kappa_points_data_list[self.current_point][1]
            for i in range(len(self.kappa_points_data_list) - 1, -1, -1):
                if self.kappa_points_data_list[i][1] > curr_ss:
                    self.current_point = i
                    break
        # down arrow key
        elif (key == 16777237):
            curr_ss = self.kappa_points_data_list[self.current_point][1]
            for i in range(len(self.kappa_points_data_list)):
                if self.kappa_points_data_list[i][1] < curr_ss:
                    self.current_point = i
                    break
        self.update_kappa_info_and_graph()

    def on_pick_kappa_points(self, event):
        """
        When mouse is clicked on a kappa point, change to focus that kappa point
        :param event:
        :return:
        """
        self.current_point = event.ind[0]
        self.update_kappa_info_and_graph()

    def toggle_exclude_include_kappa_point(self):
        """
        Toggle between including and excluding the current kappa point. If it is already excluded, then include it.
        Otherwise, exclude it.
        :return:
        """
        ss = self.kappa_points_data_list[self.current_point][1]
        dp = self.kappa_points_data_list[self.current_point][0]
        if self.kappa_points_is_included_list[(dp, ss)]:
            self.kappa_points_is_included_list[(dp, ss)] = False
        else:
            self.kappa_points_is_included_list[(dp, ss)] = True
        self.update_kappa_info_and_graph()

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
            if self.finish_scan_alignment:
                self.ccnc_sig_list = self.ccnc_sig_list_list[self.current_scan]
                self.b = self.b_list[self.current_scan]
                self.d = self.d_list[self.current_scan]
                self.c = self.c_list[self.current_scan]
                self.min_dp_asym = self.min_dp_asym_list[self.current_scan]
                self.min_dp = self.min_dp_list[self.current_scan]
                self.max_dp_asym = self.max_dp_asym_list[self.current_scan]
                self.super_saturation_rate = self.dp50_list[self.current_scan][1]
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
            self.usable_for_sigmoid_fit_list[self.current_scan] = False
            self.usable_for_kappa_cal_list[self.current_scan] = False

    def shift_data_by_one_second(self, forward=True):
        try:
            # can't shift if already pass the alignment step
            if self.finish_scan_alignment:
                return
            # invalid scan. Can't shift
            if not self.usable_for_sigmoid_fit_list[self.current_scan] or self.min_pos_SMPS_list[self.current_scan] \
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
        except:
            self.view.show_error_dialog("Can't shift processed_data. You should disable this peak!")

    def change_scan_status(self):
        # is aligning the CCNC and SMPS data of each scan
        if self.finish_scan_alignment is False:
            if self.usable_for_sigmoid_fit_list[self.current_scan] is True:
                self.usable_for_sigmoid_fit_list[self.current_scan] = False
            elif self.min_pos_CCNC_list[self.current_scan] and self.min_pos_CCNC_list[self.current_scan]:
                self.usable_for_sigmoid_fit_list[self.current_scan] = True
            else:
                self.view.show_error_dialog("You can't enable this scan! This scan is not usable!")
            self.update_view()
        # is checking sigmoid fit
        else:
            if self.usable_for_sigmoid_fit_list[self.current_scan]:
                if self.usable_for_kappa_cal_list[self.current_scan]:
                    self.usable_for_kappa_cal_list[self.current_scan] = False
                else:
                    self.usable_for_kappa_cal_list[self.current_scan] = True
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
        self.prepare_scan_data()
        self.draw_concentration_over_scan_time_graph()
        if self.finish_scan_alignment:
            if self.usable_for_sigmoid_fit_list[self.current_scan]:
                self.draw_complete_sigmoid_graph()
            else:
                self.draw_ccn_cn_ratio_over_diameter_graph()
            self.draw_all_scans_alignment_summary_graph()
            self.view.update_scan_information_after_sigmoid_fit()
        else:
            self.draw_ccn_cn_ratio_over_diameter_graph()
            self.draw_temperature_graph()
            self.view.update_scan_information()
        self.view.centralWidget().setFocus()

    def reset(self):
        i = 1
        gc.collect()
        if self.min_compare_graph:
            self.min_compare_graph.clf()
            plt.close(self.min_compare_graph)
            self.min_compare_graph = None
        if self.kappa_graph:
            self.kappa_graph.clf()
            plt.close(self.kappa_graph)
            self.kappa_graph = None
        self.current_point = None
        self.finish_scan_alignment = False
        self.cancelling_progress_bar = False
        self.flow_rate = 0.3
        self.ccn_list = None
        self.cn_list = None
        self.ccn_fixed_list = None
        self.cn_fixed_list = None
        self.g_cn_list = None
        self.g_ccn_list = None
        self.ccnc_sig_list = []
        self.ccn_normalized_list = []
        self.diameter_midpoint_list = []
        self.ccn_cn_sim_list = []
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
        self.shift_factor_list = []
        self.kappa_calculate_dict = {}
        self.alpha_pinene_dict = {}
        self.scan_position_in_data = []
        self.usable_for_kappa_cal_list = []
        self.invalid_k_points_list = []
        self.kappa_points_data_list = []
        gc.collect()

    # ---------------------progress bar-------------------------

    def move_progress_bar_forward(self, message=None, max_value=None, complete=False, value=1):
        self.view.move_progress_bar_forward(message, max_value, complete, value)
        if self.cancelling_progress_bar is True:
            self.cancelling_progress_bar = False
            raise ProgressBarInterruptException()

    def cancel_progress_bar(self):
        self.cancelling_progress_bar = True

    def run(self):
        self.parse_and_match_smps_ccnc_data()


def main():
    controller = Controller(False)
    view = View(controller)
    controller.view = view
    files = [u'C:\\Users\\KKK\\OneDrive\\Researches\\Chemics\\Examples\\AS_Calibration_SMPS.txt',
             u'C:\\Users\\KKK\\OneDrive\\Researches\\Chemics\\Examples\\CCN data 100203092813.csv']
    controller.files = files
    controller.run()
    view.show_ui()


if __name__ == '__main__':
    main()
