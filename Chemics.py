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
from HelperFunctions import *
matplotlib.style.use('ggplot')


class Controller():
    def __init__(self, mode=True):
        """"
        The controller of the program. Also works doubly as the model.
        """
        self.view = None
        self.files = None
        self.smps_txt_files = []
        self.ccnc_csv_files = []
        self.current_scan = 0
        self.number_of_scan = 0
        self.scan_start_time_list = None
        self.scan_end_time_list = None
        self.experiment_date = None
        self.time_of_each_scan = 0
        # ccnc = Cloud Condensation Nuclei Counter
        self.ccnc_data = None
        # smps = Scanning Mobility Particle Sizer
        self.smps_data = None
        self.particle_diameter_list = None
        self.dnlog_list = None
        self.ss_list = []
        self.processed_data = None
        self.raw_data = None
        self.current_point = None
        self.finish_sigmoid_fit_phase = False
        self.cancelling = False

        self.concentration = 0.3
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

        self.concentration_over_scan_time_graph = None
        self.complete_dry_diameter_graph = None
        self.min_compare_grpah = None
        self.kappa_graph = None
        self.temperature_graph = None
        self.adjusted_graph_list = []
        self.dry_diameter_graph_list = []
        self.min_pos_SMPS_list = []
        self.min_pos_CCNC_list = []
        self.temperature_graph_list = []

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
        self.kappa_exclude_list = []
        self.kappa_points = []
        self.klines = None
        self.max_kappa = None
        self.min_kappa = None
        self.is_full_kappa_graph = False
        self.kappa_calculate_dict = {}
        self.alpha_pinene_dict = {}
        self.scan_position_in_data = []
        self.usable_for_kappa_cal_list = []
        self.clone = None


    # -------------General processed_data processing------------------

    def get_concentration(self):
        """
        """
        # TODO: what is the concentration here? ask the professor
        self.concentration = self.view.get_concentration()

    def get_raw_data(self):
        """
        Get the smps and ccnc processed_data from the input files.
        """
        self.move_progress_bar_forward("Reading data files...")
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
        self.ccnc_data = process_csv_files(self.ccnc_csv_files)
        self.smps_data = process_text_files(self.smps_txt_files)
        self.experiment_date = self.ccnc_data[0]
        self.scan_start_time_list = self.smps_data[0]

        # Get the scan time (in seconds) of each scan
        for i in range(len(self.smps_data)):
            if ''.join(self.smps_data[i][0].split()).lower() == "scanuptime(s)":
                scan_up_time = self.smps_data[i][1]
                scan_down_time = self.smps_data[i+1][1]
                break
        self.time_of_each_scan = int(scan_up_time) + int(scan_down_time)
        self.scan_end_time_list = []
        for i in range(len(self.scan_start_time_list)):
            self.scan_end_time_list.append(datetime.strptime(self.scan_start_time_list[i], "%H:%M:%S") + timedelta(
                seconds=self.time_of_each_scan))

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

    def get_dnlog_data(self):
        # TODO: what is the dnLog
        self.move_progress_bar_forward("Acquiring dnlog data...")
        for k in range(1, len(self.smps_data)):
            if re.search('[aParam-zA-Z]', self.smps_data[k][0]):
                dnlogPos = k
                break
        startTime = ['dp'] + self.scan_start_time_list
        endTime = ['dp'] + self.scan_end_time_list
        # TODO: what is dndlogdp???????
        self.dnlog_list = [startTime, endTime]
        for k in range(1, dnlogPos):
            self.dnlog_list.append(self.smps_data[k][:len(startTime)])

    def get_smps_and_ccnc_data(self):
        self.move_progress_bar_forward("Acquiring SMPS and CCNC data...")
        ccnc_data = self.ccnc_data[1]
        smps_data = self.smps_data
        startTime = self.scan_start_time_list
        endTime = self.scan_end_time_list
        width = self.number_of_scan

        # get the count position
        for k in range(3, len(self.smps_data)):
            if re.search('[aParam-zA-Z]', self.smps_data[k][0]):
                for j in range(k, len(self.smps_data)):
                    if not re.search('[aParam-zA-Z]', self.smps_data[j][0]):
                        countPos = j
                        break
                break

        # Get the count for each scan smps
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

        # Get the CCNC processed_data
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

        # Get the first time stamp position in the ccnc file
        while (True):
            timeStamp = startTime[0]
            if ccnc_data[0][0] > timeStamp:
                startTime = startTime[1:]
                endTime = endTime[1:]
            else:
                break
        for k in range(0, len(ccnc_data)-1):
            if ccnc_data[k][0] == timeStamp and ccnc_data[k+1][0] != timeStamp:
                break

        i = 0
        t = 0
        isPeakTimeStamp = True
        smpsCcnList = []
        timeStamp = datetime.strptime(startTime[0], "%H:%M:%S")
        self.scan_position_in_data = []
        while (True):
            sizeSum = 0
            countSum = 0
            previousTimeStamp = datetime.strptime(startTime[i], "%H:%M:%S").time()
            maxSS = 0
            # Get all processed_data within the time frame
            self.scan_position_in_data.append(len(smpsCcnList))
            for t in range(self.time_of_each_scan):
                aLine = [smpsList[t][1]] + [timeStamp.time()] + [float(smpsList[t][0])] + [
                    float(smpsList[t][i + 2])]
                sizeSum = 0
                countSum = 0
                aLine.append(float(ccnc_data[k + t][-3]))
                if ccnc_data[k + t][1] > maxSS:
                    maxSS = ccnc_data[k + t][1]
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
            k += self.time_of_each_scan

            # Update the lists
            self.ss_list.append(maxSS)

            # If reach the end of smps, collect the rest of ccnc for alignment
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
                startTime[i+1] = (datetime.strptime(startTime[i + 1], "%H:%M:%S") + timedelta(seconds = abs(timeGap))).strftime( "%I:%M:%S")
                self.scan_start_time_list = startTime
                timeGamp = 0
            for t in range(timeGap):
                aLine = [None] + [timeStamp.time()] + [None] + [None]
                sizeSum = 0
                countSum = 0
                aLine.append(float(ccnc_data[k + t][-3]))
                if ccnc_data[k + t][1] > maxSS:
                    maxSS = ccnc_data[k + t][1]
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

        title = ["scan time"] + ["real time"] + ["dp"] + ["SMPS Count"] + ["CCNC Count"] + ["Ave Size"] + ["T1"] + ["T2"] + ["T3"]
        self.raw_data = pandas.DataFrame(smpsCcnList, columns=title)

    def match_smps_ccnc_data(self):
        self.move_progress_bar_forward("Matching SMPS and CCNC data....")
        self.usable_for_kappa_cal_list = []
        self.min_pos_CCNC_list = []
        self.min_pos_SMPS_list = []
        additional_data_count = int(self.time_of_each_scan / 2)
        start_time = 0
        end_time = self.time_of_each_scan
        current_scan = 0
        new_data = [self.raw_data.columns.values.tolist()]
        min_dist = self.time_of_each_scan / 10
        shift_factor = 0
        self.shift_factor_list = []
        self.number_of_scan = 0
        while True:
            a_scan = numpy.asarray(self.raw_data[start_time:end_time]["SMPS Count"])
            min_pos_smps = get_min_index(a_scan)
            if min_pos_smps == -1:
                self.min_pos_SMPS_list.append(None)
                self.min_pos_CCNC_list.append(None)
                self.usable_for_kappa_cal_list.append(False)
                min_pos_smps = 0
                min_pos_ccnc = 0
            else:
                self.min_pos_SMPS_list.append(min_pos_smps + start_time)
                a_scan = numpy.asarray(self.raw_data[start_time + shift_factor:end_time + additional_data_count]["CCNC Count"])
                min_pos_ccnc = get_min_index(a_scan)
                if min_pos_ccnc == -1:
                    self.min_pos_SMPS_list[-1] = None
                    self.min_pos_CCNC_list.append(None)
                    self.usable_for_kappa_cal_list.append(False)
                    min_pos_smps = 0
                    min_pos_ccnc = 0
                else:
                    self.usable_for_kappa_cal_list.append(True)
                    shift_factor = shift_factor + min_pos_ccnc - min_pos_smps
                    if current_scan == 0:
                        self.min_pos_CCNC_list.append(min_pos_smps)
                    else:
                        self.min_pos_CCNC_list.append(min_pos_ccnc + start_time + (min_pos_ccnc - min_pos_smps))
            self.shift_factor_list.append(shift_factor)

            for i in range(self.time_of_each_scan):
                scan_time = i + 1
                real_time = self.raw_data.iat[start_time + i, 1]
                dp = self.raw_data.iat[start_time + i, 2]
                smpsCount = self.raw_data.iat[start_time + i, 3]
                ccncCount = self.raw_data.iat[start_time + i + shift_factor, 4]
                aveSize = self.raw_data.iat[start_time + i + shift_factor, 5]
                t1 = self.raw_data.iat[start_time + i + shift_factor, 6]
                t2 = self.raw_data.iat[start_time + i + shift_factor, 7]
                t3 = self.raw_data.iat[start_time + i + shift_factor, 8]
                new_data.append([scan_time, real_time, dp, smpsCount, ccncCount, aveSize,t1,t2,t3])

            current_scan += 1
            self.number_of_scan += 1
            if end_time + self.time_of_each_scan + shift_factor >= len(self.raw_data) or current_scan >= len(self.scan_position_in_data):
                break

            start_time = self.scan_position_in_data[current_scan]
            end_time = start_time + self.time_of_each_scan
        headers = new_data.pop(0)
        self.processed_data = pandas.DataFrame(new_data, columns=headers)
        self.processed_data = convert_date(self.processed_data)

    def create_all_scan_alignment_summary_graph(self):
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
        # Plot only valid points
        for i in range(len(self.min_pos_CCNC_list)):
            if self.min_pos_CCNC_list[i] and self.min_pos_SMPS_list[i]:
                tempSMPSPeakList.append(self.min_pos_SMPS_list[i])
                tempCCNPeakCList.append(self.min_pos_CCNC_list[i])
            else:
                # Make up for the null values
                if len(tempSMPSPeakList) >0:
                    tempSMPSPeakList.append(tempSMPSPeakList[-1] + self.time_of_each_scan)
                    tempCCNPeakCList.append(tempCCNPeakCList[-1] + self.time_of_each_scan)
                else:
                    tempSMPSPeakList.append(10)
                    tempCCNPeakCList.append(10)

        x = numpy.asarray(tempSMPSPeakList)
        y = numpy.asarray(tempCCNPeakCList)

        result = scipy.stats.linregress(x, y)
        slope = result[0]
        yIntercept = result[1]

        # Recalculate the position of the smps
        plt.plot(x, x * slope + yIntercept, linewidth = 4, color = '#43A047', label = "Regression line")
        plt.plot(x, y, "o", ms = 10, color = "#43A047",picker=5, mew = 0, label = "Minimum")
        textToShow = str(slope) + "* x" + " + " + str(yIntercept)
        # plt.text(x[4], y[3], textToShow, color="#81C784")
        self.current_point, = plt.plot(x[0],y[0],'o', color = "#81C784", ms = 12, mew = 0)
        plt.xlabel("SMPS minumum point")
        plt.ylabel("CCNC minimum point")
        handles, labels = plt.gca().get_legend_handles_labels()
        legend = plt.legend(handles, labels,loc="upper right", bbox_to_anchor=(1,0.7))
        legend.get_frame().set_facecolor('#9E9E9E')

        if not self.min_compare_grpah:
            self.min_compare_grpah = plt.gcf()
        else:
            plt.close(self.min_compare_grpah)
            self.min_compare_grpah = plt.gcf()

    def on_pick(self, event):
        """
        When aParam dot on the graph is selected
        """
        peak = event.ind[0]
        if peak != self.current_scan:
            self.switch_to_scan(peak)

    # -------------Single peak processing ------------------

    def prepare_scan_data(self):
        try:
            dNdLogDpList = []
            startPoint = self.current_scan * self.time_of_each_scan
            endPoint = (self.current_scan + 1) * self.time_of_each_scan
            self.particle_diameter_list = list(self.processed_data[startPoint:endPoint]['dp'])
            self.ccn_list = list(self.processed_data[startPoint:endPoint]['CCNC Count'])
            self.cn_list = list(self.processed_data[startPoint:endPoint]['SMPS Count'])
            self.temp1 = list(self.processed_data[startPoint:endPoint]['T1'])
            self.temp2 = list(self.processed_data[startPoint:endPoint]['T2'])
            self.temp3 = list(self.processed_data[startPoint:endPoint]['T3'])

            # modify this to get aParam decently good number
            self.cn_list = [x * 0.2 for x in self.cn_list]
            self.diameter_midpoint_list = []
            self.ccn_cn_ratio_list = []
            for i in range(len(self.ccn_list)):
                self.ccn_cn_ratio_list.append(self.ccn_list[i] / self.cn_list[i])
            self.drop_size_list = list(self.processed_data[startPoint:endPoint]['Ave Size'])
            for i in range(2, len(self.dnlog_list)):
                self.diameter_midpoint_list.append(self.dnlog_list[i][0])
                dNdLogDpList.append(self.dnlog_list[i][self.current_scan + 1])
            self.ccn_normalized_list = normalize_list(dNdLogDpList)
            self.usable_for_kappa_cal_list[self.current_scan] = True
        except:
            self.min_pos_CCNC_list[self.current_scan] = None
            self.min_pos_SMPS_list[self.current_scan] = None
            self.usable_for_kappa_cal_list[self.current_scan] = False

    def init_correct_charges(self):
        """
        Initiate the correct charge procedure
        """
        # If the peak doesn't have good alignment, do nothing
        self.cn_fixed_list = self.cn_list[:]
        self.ccn_fixed_list = self.ccn_list[:]
        self.g_ccn_list = self.ccn_fixed_list[:]
        self.g_cn_list = self.cn_fixed_list[:]

    def correct_charges(self):
        try:
            self.move_progress_bar_forward()
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
                QtGui.qApp.processEvents()
                aDList = [0]
                for k in range(1, 4):
                    c = cal_cc(i * 10 ** -9, lambdaAir)
                    dp = 10 ** 9 * find_dp(i * 10 ** -9 / c, lambdaAir, k)
                    aDList.append(dp)
                chargeList.append(aDList)

            # second part of correct charges
            self.cn_fixed_list = self.cn_list[:]
            self.ccn_fixed_list = self.ccn_list[:]
            maxUpperBinBound = (self.particle_diameter_list[-1] + self.particle_diameter_list[-2]) / 2
            lenDpList = len(self.particle_diameter_list)

            for i in range(lenDpList):
                QtGui.qApp.processEvents()
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
                                    self.ccn_fixed_list[j] = self.ccn_fixed_list[j] + moveDoubletCounts * self.g_ccn_list[j] / \
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
                                self.ccn_fixed_list[j] = self.ccn_fixed_list[j] + moveTripletCounts * self.ccn_list[j] / \
                                                                              self.cn_list[j]
                            else:
                                self.ccn_fixed_list[j] = self.ccn_fixed_list[j] + moveTripletCounts
                            break
                        j -= 1

            for i in range(len(self.ccn_fixed_list)):
                QtGui.qApp.processEvents()
                if self.ccn_fixed_list[i] / self.cn_fixed_list[i] < -0.01:
                    self.ccn_fixed_list[i] = 0

            self.g_ccn_list = self.ccn_fixed_list[:]
            self.g_cn_list = self.cn_fixed_list[:]
            self.usable_for_kappa_cal_list[self.current_scan] = True
        except:
            raise ScanDataError()

    def get_constants(self, min_dp = None, min_dp_asym = None, max_dp_asym = None):
        """
        Acquire the necessary constants from the processed_data to perform sigmodal fit
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
                checkLength = self.time_of_each_scan / 20
                for i in range(checkLength, len(self.ccnc_sig_list)-checkLength):
                    isMid = False
                    if self.ccnc_sig_list[i] > 0.5:
                        isMid = True
                        # check if the previous 5 numbers are smaller and the next 5 numbers are bigger
                        for j in range(1,checkLength):
                            if self.ccnc_sig_list[i+j] < self.ccnc_sig_list[i]:
                                isMid = False
                                break
                    if isMid:
                        mpIndex = i
                        break

                minDpAsymPos = 0
                currMax = 0
                # Get minDP
                for i in range(mpIndex,1, -1):
                    if self.ccnc_sig_list[i] < 0.1:
                        self.min_dp = self.particle_diameter_list[i]
                        break

                # Get minDpAsym
                for i in range(mpIndex, mpIndex + self.time_of_each_scan / 10):
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
            self.ccn_cn_sim_list.append(0)
            for i in range(1, len(self.particle_diameter_list)):
                if self.min_dp < self.particle_diameter_list[i] < self.max_dp:
                    n = self.b / (1 + (self.particle_diameter_list[i] / self.d) ** self.c)
                    self.ccn_cn_sim_list.append(n)
                else:
                    self.ccn_cn_sim_list.append(self.ccn_cn_sim_list[i - 1])
            self.usable_for_kappa_cal_list[self.current_scan] = True
        except:
            raise ScanDataError()

    def fit_sigmoid_line_to_run(self):
        """
        fit the sigmoid line to the the data points.
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

            for i in range(1, len(self.particle_diameter_list)):
                if self.min_dp < self.particle_diameter_list[i] < self.max_dp:
                    n = self.b / (1 + (self.particle_diameter_list[i] / self.d) ** self.c)
                    self.ccn_cn_sim_list.append(n)
                else:
                    self.ccn_cn_sim_list.append(self.ccn_cn_sim_list[i - 1])
            self.usable_for_kappa_cal_list[self.current_scan] = True
        except:
            raise ScanDataError()

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
            x = range(self.time_of_each_scan)
            plt.plot(x, self.cn_list, linewidth = 4, color = '#EF5350', label = "CN")
            plt.plot(x, self.ccn_list, linewidth = 4, color = '#2196F3',label = "CCN")
            handles, labels = plt.gca().get_legend_handles_labels()
            legend = plt.legend(handles, labels, loc="upper left", bbox_to_anchor=(0, 0.9))
            legend.get_frame().set_facecolor('#9E9E9E')

            plt.xlabel("Scan time(s)")
            plt.ylabel("Concentration(cm3)")
        except:
            figure = plt.figure(facecolor=settings.graphBackgroundColor)
        finally:
            self.adjusted_graph_list.append(plt.gcf())

    def create_ccn_cn_ratio_over_diameter_graph(self, new_figure=one):
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
            plt.axhline(1, color='0.7', linewidth=2, linestyle = 'dashed')
            plt.gca().tick_params(axis='x', color='1', which='both', labelcolor="0.6")
            plt.gca().tick_params(axis='y', color='1', which='both', labelcolor="0.6")
            plt.gca().yaxis.label.set_color('0.6')
            plt.gca().xaxis.label.set_color('0.6')
            plt.plot(self.particle_diameter_list, self.ccn_cn_ratio_list, 'o', color="#2196F3", mew=0.5, mec ="#0D47A1", ms = 9, label ="CCN/CN")
            yLim = min(2, max(self.ccn_cn_ratio_list)) + 0.2
            plt.gca().axes.set_ylim([-0.1, yLim])
            handles, labels = plt.gca().get_legend_handles_labels()
            legend = plt.legend(handles, labels,loc="upper right", bbox_to_anchor=(1,0.9))
            legend.get_frame().set_facecolor('#9E9E9E')
            plt.xlabel("Diameter (nm)")
            plt.ylabel("CCN/CN")
        except:
            figure = plt.figure(facecolor=settings.graphBackgroundColor)
        finally:
            self.dry_diameter_graph_list.append(plt.gcf())

    def create_temperature_graph(self, new_figure = None):
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

            x = range(self.time_of_each_scan)
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

    def create_complete_sigmoid_graph(self, new_figure = None):
        """
        Make complete graph of the dry diameter after optimization and sigmodal fit
        """
        try:
            if not self.min_pos_CCNC_list[self.current_scan] or not self.min_pos_SMPS_list[self.current_scan]:
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
            plt.axhline(1, color='0.7', linewidth=2, linestyle = "--")
            plt.gca().tick_params(axis='x', color='1', which='both', labelcolor="0.6")
            plt.gca().tick_params(axis='y', color='1', which='both', labelcolor="0.6")
            plt.gca().yaxis.label.set_color('0.6')
            plt.gca().xaxis.label.set_color('0.6')
            yLim = min(2, max(self.ccnc_sig_list)) + 0.2
            plt.gca().axes.set_ylim([-0.1, yLim])
            plt.plot(self.diameter_midpoint_list, self.ccn_normalized_list, linewidth=4, color='#43A047', label="dN/dLogDP")
            if self.usable_for_kappa_cal_list[self.current_scan]:
                plt.plot(self.particle_diameter_list, self.ccn_cn_sim_list, linewidth=5, color='#EF5350', label="Sigmodal Fit")
            plt.plot(self.particle_diameter_list, self.ccn_cn_ratio_list, 'o', color="#2196F3", mew=0.5, mec ="#1976D2", ms = 9, label ="CCN/CN")
            plt.plot(self.particle_diameter_list, self.ccnc_sig_list, 'o', color="#1565C0", mew=0.5, mec ="#0D47A1", ms = 9, label ="CCN/CN (Corrected)")
            plt.xlabel("Dry diameter(nm)")
            plt.ylabel("CCN/CN ratio and Normalized dN/dLogDP")
            handles, labels = plt.gca().get_legend_handles_labels()
            legend = plt.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.7, 1.1))
            legend.get_frame().set_facecolor('#9E9E9E')
        except:
            figure = plt.figure(facecolor=settings.graphBackgroundColor)
        finally:
            self.dry_diameter_graph_list.append(plt.gcf())

    def create_kappa_graph(self):
        """
        Produce the kappa graph, may be in full or only around the points
        """
        # Read in the kappa lines csv files.
        if self.klines is None:
            self.klines = pandas.read_csv("klines.csv", header=1)
        # Acquire the kappa points
        kappaList = []
        stdKappaList = []
        self.kappa_points = []

        for key in self.alpha_pinene_dict.keys():
            kappaList.append(self.alpha_pinene_dict[key][2])
            stdKappaList.append(self.alpha_pinene_dict[key][3])
            if not self.is_full_kappa_graph:
                self.kappa_points.append((self.alpha_pinene_dict[key][0], key))
            else:
                for i in self.alpha_pinene_dict[key][-1]:
                    self.kappa_points.append((i, key))

        # Determine maximum and minimum kappas
        temp_kappa_list = []
        for i in range(len(kappaList)):
            temp_kappa_list.append(kappaList[i] + stdKappaList[i])
        self.max_kappa = max(temp_kappa_list)
        temp_kappa_list = []
        for i in range(len(kappaList)):
            temp_kappa_list.append(kappaList[i] - stdKappaList[i])
        self.min_kappa = min(temp_kappa_list)
        header = self.klines.columns
        diaList = self.klines[header[1]]

        fullKXList = []
        fullKYList = []
        kpXList = []
        kpYList = []
        excludedXList = []
        excludedYList = []

        for i in range(len(self.kappa_points)):
            fullKXList.append(self.kappa_points[i][0])
            fullKYList.append(self.kappa_points[i][1])
        # Process which point is included and which point is excluded
        # Include all excluded points to exclude list
        for i in self.kappa_exclude_list:
            excludedXList.append(self.kappa_points[i][0])
            excludedYList.append(self.kappa_points[i][1])
        # Exclude all excluded points from kpx/y list
        for i in range(len(self.kappa_points)):
            if i not in self.kappa_exclude_list:
                kpXList.append(self.kappa_points[i][0])
                kpYList.append(self.kappa_points[i][1])
        # convert list to array
        excludedXList = numpy.asarray(excludedXList)
        excludedYList = numpy.asarray(excludedYList)
        kpXList = numpy.asarray(kpXList)
        kpYList = numpy.asarray(kpYList)

        # Get all kappa points - without average
        # Prepare the figure
        if self.kappa_graph:
            figure = plt.figure(self.kappa_graph.number)
            figure.clf()
        else:
            figure = plt.figure(facecolor=settings.graphBackgroundColor)
        plt.axes(frameon=False)
        plt.grid(color='0.5')
        plt.axhline(0.1, color='0.6', linewidth=4)
        plt.axvline(10, color='0.6', linewidth=4)
        plt.gca().tick_params(axis='x', color='1', which='both', labelcolor="0.6")
        plt.gca().tick_params(axis='y', color='1', which='both', labelcolor="0.6")
        plt.gca().yaxis.label.set_color('0.6')
        plt.gca().xaxis.label.set_color('0.6')
        plt.ylim([0.1, 1.5])
        plt.xlim([10, 200])
        plt.grid(True, which='both', color="0.85")
        plt.xlabel("Dry diameter(nm)")
        plt.ylabel("Super Saturation(%)")
        figure.canvas.mpl_connect('pick_event', self.on_kappa_pick)

        # Graph the klines
        # Full kappa lines
        #if self.is_full_kappa_graph:
         #   for i in range(2, len(header)):
         #       y = self.klines[header[i]]
        #        plt.loglog(diaList, y, label=str(header[i]), linewidth=4)
        # Draw only the portion around the kappa
       # else:
        i = 2
        kappa = 1
        step = 0.1
        kappaStartPos = 2
        kappaEndPos = len(header)
        while True:
            if self.max_kappa > kappa:
                kappaStartPos = max(2, i - 3)
                break
            i += 1
            kappa -= step
            if kappa == step:
                step /= 10
            if i >= len(header):
                kappaStartPos = len(header)
                break
        i = 2
        kappa = 1
        step = 0.1
        while True:
            if self.min_kappa > kappa:
                kappaEndPos = min(i + 3, len(header))
                break
            i += 1
            kappa -= step
            if kappa == step:
                step /= 10
            if i >= len(header):
                kappaEndPos = len(header)
                break
        for i in range(kappaStartPos, kappaEndPos):
            y = self.klines[header[i]]
            plt.loglog(diaList, y, label=str(header[i]), linewidth=4)

        # Graph all the kappa points
        plt.plot(fullKXList, fullKYList, 'o', picker=5, mew=0.5, ms=12, alpha = 0)
        # Graph the kappa points
        plt.plot(kpXList, kpYList, 'o', color="#81C784", mew=0.5, mec="#81C784", ms=12,
                 label=" Kappa Points")
        # Graph the excluded points
        if len(excludedXList) >0 and len(excludedYList) > 0:
            plt.plot(excludedXList, excludedYList, 'x', mec="#81C784", color="#81C784",mew=4,ms = 12, label="excluded kappa")
        handles, labels = plt.gca().get_legend_handles_labels()
        legend = plt.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.1, 1))
        legend.get_frame().set_facecolor('#9E9E9E')
        self.kappa_graph = plt.gcf()
        self.view.update_kappa_graph()

    def on_kappa_pick(self, event):
        """
        When a kappa point is clicked
        """
        kappaPoint = event.ind[0]
        excluded = False
        # if already in excluded list, include the points
        for i in range(len(self.kappa_exclude_list)):
            if self.kappa_exclude_list[i] == kappaPoint:
                self.kappa_exclude_list = self .kappa_exclude_list[:i] + self.kappa_exclude_list[i + 1:]
                excluded = True
                break
        # else, exclude the point
        if not excluded:
            self.kappa_exclude_list.append(kappaPoint)
        # Remake the kappa graph
        self.create_kappa_graph()

    def correct_charges_and_fit_sigmoid_one_scan(self):
        """
        optimize a single run
        """
        try:
            # If a peak is invalid, then do nothing
            if not self.min_pos_CCNC_list[self.current_scan] or not self.min_pos_SMPS_list[self.current_scan]:
                self.prepare_scan_data()
                raise ScanDataError()
            else:
                self.ccnc_sig_list = []
                self.prepare_scan_data()
                self.move_progress_bar_forward()
                remove_small_ccn(self.ccn_list, self.min_ccn)
                self.move_progress_bar_forward()
                self.init_correct_charges()
                for i in range(5):
                    self.correct_charges()
                for i in range(len(self.ccn_fixed_list)):
                    self.ccnc_sig_list.append(self.ccn_fixed_list[i] / self.cn_fixed_list[i])
                self.move_progress_bar_forward()
                self.get_constants()
                self.move_progress_bar_forward()
                self.fit_sigmoid_line_to_run()
                self.move_progress_bar_forward()
                plt.ioff()
                self.create_complete_sigmoid_graph()
                self.move_progress_bar_forward()
        except ScanDataError:
            # Disable the peak
            self.usable_for_kappa_cal_list[self.current_scan] = False
            # Add empty processed_data to the list
            self.b_list.append(0)
            self.d_list.append(0)
            self.c_list.append(0)
            self.dp50_list.append((0, 0))
            self.dp50_wet_list.append(0)
            self.dp50_plus_20_list.append(0)
            self.min_dp_list.append(0)
            self.min_dp_asym_list.append(0)
            self.max_dp_asym_list.append(0)
            self.create_complete_sigmoid_graph()
        else:
            # Store processed_data
            self.min_dp_list.append(self.min_dp)
            self.min_dp_asym_list.append(self.min_dp_asym)
            self.max_dp_asym_list.append(self.max_dp_asym)
            self.b_list.append(self.b)
            self.d_list.append(self.d)
            self.c_list.append(self.c)
            self.dp50_wet_list.append(self.dp50_wet)
            self.dp50_plus_20_list.append(self.dp50_plus_20)
            self.dp50_list.append((self.d, self.ss_list[self.current_scan]))

    def parse_and_match_smps_ccnc_data(self):
        try:
            self.reset()
            self.view.reset()
            self.move_progress_bar_forward(max_value=5)
            self.get_concentration()
            self.get_raw_data()
            self.get_dnlog_data()
            self.get_smps_and_ccnc_data()
            self.match_smps_ccnc_data()
            self.move_progress_bar_forward(complete=True)
            self.move_progress_bar_forward(max_value=self.number_of_scan)
            for i in range(0, self.number_of_scan):
                self.move_progress_bar_forward("Processing scan number " + str(i + 1), value=0)
                plt.ioff()
                self.current_scan = i
                self.prepare_scan_data()
                self.create_concentration_over_scan_time_graph()
                self.create_ccn_cn_ratio_over_diameter_graph()
                self.create_temperature_graph()
            self.move_progress_bar_forward(complete=True)
            self.current_scan = -1
            self.view.update_general_information()
            self.create_all_scan_alignment_summary_graph()
            self.switch_to_scan(0)
        except FileNotFoundError:
            self.view.show_error_dialog("Can't find SMPS or CCNC processed_data in files: " + str(self.files))
            raise DataPreparationError()
        except FileProcessingError:
            self.view.show_error_dialog("Can't process the SMPS or CCNC files!")
            raise DataPreparationError()
        except dnlogDataError:
            self.view.show_error_dialog("Can't process dnlog processed_data from the SMPS file!")
            raise DataPreparationError()
        except DataError:
            self.view.show_error_dialog("Can't process processed_data from SMPS or CCNC raw processed_data!")
            raise DataPreparationError()
        except DataMatchingError:
            self.view.show_error_dialog("Can't match SMPS and CCNC processed_data!")
            raise DataPreparationError()
        except InterruptError:
            self.reset()
            self.view.reset()
            self.view.show_error_dialog("The preparation process is cancelled")
            raise DataPreparationError()

    def correct_charges_and_fit_sigmoid_all_scans(self):

        for aFigure in self.temperature_graph_list:
            aFigure.clf()
            plt.close(aFigure)
        self.temperature_graph_list = []
        for aFigure in self.dry_diameter_graph_list:
            aFigure.clf()
            plt.close(aFigure)
        self.dry_diameter_graph_list = []

        try:
            self.move_progress_bar_forward(max_value=self.number_of_scan * 12)
            self.finish_sigmoid_fit_phase = True
            for i in range(0, self.number_of_scan):
                self.current_scan = i
                self.move_progress_bar_forward("Correcting charges and fitting sigmoid for scan number " + str(i + 1), value=0)
                self.correct_charges_and_fit_sigmoid_one_scan()
            self.current_scan = -1
            self.move_progress_bar_forward(complete=True)
            self.switch_to_scan(0)
        except InterruptError:
            self.view.show_error_dialog("The sigmoid fitting process is cancelled!")
            self.finish_sigmoid_fit_phase = False
            # TODO: make a reroll mechanism to allow going back to previous step when there is an interruption
            self.parse_and_match_smps_ccnc_data()


    #----------------Kappa Calculation--------------------------

    def setConstantForCalKappa(self, sigma=None, temp=None, dd1=None, iKappa=None, dd2=None, iKappa2=None, solu=None):
        """
        Update the necessary constants to calculate Kappa
        """
        if sigma:
            self.sigma = sigma
        if temp:
            self.temp = temp
        if dd:
            self.dd = dd
        if iKappa:
            self.iKappa = iKappa
        if dd2:
            self.dd2 = dd2
        if iKapp2:
            self.iKappa2 = iKappa2
        if solu:
            self.solubility = solu

    def calKappa(self):
        """
        Calculate the kappa values - producing both raw processed_data kappa and graph processed_data kappa
        """
        for i in range(len(self.dp50_list)):
            self.usable_for_kappa_cal_list.append(True)
        self.move_progress_bar_forward("Calculating Kappa...", max_value=len(self.dp50_list) + 3)
        self.move_progress_bar_forward("Reading in lookup table")
        if self.kappa_excel is None:
            self.kappa_excel = pandas.read_csv("kCal.csv", header=None)
        lookup = self.kappa_excel
        self.move_progress_bar_forward("Calculating basic consts...")
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
                self.move_progress_bar_forward("Calculating sc...")
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
            self.move_progress_bar_forward("Processing peak" + str(i + 1))
            self.trueSC = (max(sList1[i:]) - 1) * 100
            self.anaKappa = (4 * self.aParam ** 3) / (27 * (dp50 * 0.000000001) ** 3 * log(ss / 100 + 1) ** 2)
            kDevi = (self.appKappa - self.anaKappa) / self.appKappa * 100
            if ss in self.kappa_calculate_dict.keys():
                self.kappa_calculate_dict[ss].append([dp50, self.appKappa, self.anaKappa, kDevi, self.trueSC])
            else:
                self.kappa_calculate_dict[ss] = ([[dp50, self.appKappa, self.anaKappa, kDevi, self.trueSC]])

        for aKey in self.kappa_calculate_dict.keys():
            aSSList = self.kappa_calculate_dict[aKey]
            dp50List = []
            appKappaList = []
            anaKappaList = []
            meanDevList = []
            for aSS in aSSList:
                dp50List.append(aSS[0])
                appKappaList.append(aSS[1])
                anaKappaList.append(aSS[2])
                meanDevList.append(aSS[3])
            meanDp = average(dp50List)
            stdDp = numpy.std(dp50List)
            meanApp = average(appKappaList)
            stdApp = numpy.std(appKappaList)
            meanAna = average(anaKappaList)
            stdAna = numpy.std(anaKappaList)
            meanDev = average(meanDevList)
            devMean = (meanApp - meanAna) / meanApp * 100
            self.alpha_pinene_dict[aKey] = (meanDp, stdDp, meanApp, stdApp, meanAna, stdAna, meanDev, devMean, dp50List)
        self.move_progress_bar_forward(complete=True)

    #--------------- UI interaction - Data interaction functions--------------------

    def updateGUI(self):
        """
        Update the GUI with the newest processed_data of the current peak
        """
        self.concentration_over_scan_time_graph = self.adjusted_graph_list[self.current_scan]
        self.complete_dry_diameter_graph = self.dry_diameter_graph_list[self.current_scan]
        self.view.update_dp_dnlog_figures(self.concentration_over_scan_time_graph, self.complete_dry_diameter_graph)
        if self.finish_sigmoid_fit_phase:
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
            self.current_point.set_xdata(numpy.asarray(self.min_pos_SMPS_list)[self.current_scan])
            self.current_point.set_ydata(numpy.asarray(self.min_pos_CCNC_list)[self.current_scan])
            self.view.update_temp_and_min_figure(self.min_compare_grpah)
            self.view.update_information_after_sig_fit()

        else:
            self.temperature_graph = self.temperature_graph_list[self.current_scan]
            self.super_saturation_rate = self.ss_list[self.current_scan]
            self.view.update_temp_and_min_figure(self.temperature_graph)
            self.view.update_information_of_run()

    def shift_ccnc_data_by_one_second(self, forward=True):
        try:
            if not self.usable_for_kappa_cal_list[self.current_scan] or self.min_pos_SMPS_list[self.current_scan] is None or \
                            self.min_pos_CCNC_list[self.current_scan] is None:
                return
            shift_factor = self.shift_factor_list[self.current_scan]
            if forward is True:
                shift_factor -= 1
            else:
                shift_factor += 1
            self.shift_factor_list[self.current_scan] = shift_factor
            # Get the new processed_data list from raw processed_data
            start_time = self.scan_position_in_data[self.current_scan] + shift_factor
            end_time = start_time + self.time_of_each_scan
            new_ccn = list(self.raw_data.iloc[start_time:end_time, 4])
            new_ave_size = list(self.raw_data.iloc[start_time:end_time, 5])
            new_t1 = list(self.raw_data.iloc[start_time:end_time, 6])
            new_t2 = list(self.raw_data.iloc[start_time:end_time, 7])
            new_t3 = list(self.raw_data.iloc[start_time:end_time, 8])

            # Get the location of the replacement in the processed_data
            data_start_time = self.time_of_each_scan * self.current_scan
            data_end_time = data_start_time + self.time_of_each_scan

            # Replace the processed_data in the processed_data
            for i in range(self.time_of_each_scan):
                self.processed_data.iat[data_start_time + i, 4] = new_ccn[i]
                self.processed_data.iat[data_start_time + i, 5] = new_ave_size[i]
                self.processed_data.iat[data_start_time + i, 6] = new_t1[i]
                self.processed_data.iat[data_start_time + i, 7] = new_t2[i]
                self.processed_data.iat[data_start_time + i, 8] = new_t3[i]

            self.prepare_scan_data()
            figure = self.concentration_over_scan_time_graph
            self.create_concentration_over_scan_time_graph(new_figure=figure)
            figure = self.complete_dry_diameter_graph
            self.create_ccn_cn_ratio_over_diameter_graph(figure)
            figure = self.temperature_graph
            self.create_temperature_graph(figure)
            self.view.updateDpdnlogFigures(self.concentration_over_scan_time_graph, self.complete_dry_diameter_graph)
            self.view.update_temp_and_min_figure(self.temperature_graph)
        except:
            self.view.show_error_dialog("Can't shift processed_data. You should disable this peak!")

    def disable_scan(self):
        if not self.finish_sigmoid_fit_phase:
            self.min_pos_SMPS_list[self.current_scan] = None
            self.min_pos_CCNC_list[self.current_scan] = None
            self.usable_for_kappa_cal_list[self.current_scan] = False
        else:
            self.usable_for_kappa_cal_list[self.current_scan] = False
        self.updateGUI()

    def switch_to_scan(self, peak):
        """
        Switch to a new peak
        :param peak: the new peak
        """
        if peak != self.current_scan:
            self.current_scan = peak
            self.updateGUI()

    def refitting_sigmoid_line(self, min_dry_diameter, min_dry_diameter_asymptote, max_dry_diameter_asymptote):
        if not self.finish_sigmoid_fit_phase:
            return
        self.move_progress_bar_forward("Updating sigmoid line fitting of current scan..." + str(self.current_scan + 1), max_value=6)
        self.usable_for_kappa_cal_list[self.current_scan] = True
        self.ccnc_sig_list = []
        try:
            self.prepare_scan_data()
            self.move_progress_bar_forward()
            remove_small_ccn(self.ccn_list, self.min_ccn)
            self.move_progress_bar_forward()
            self.init_correct_charges()
            for i in range(5):
                self.correct_charges()
            for i in range(len(self.ccn_fixed_list)):
                self.ccnc_sig_list.append(self.ccn_fixed_list[i] / self.cn_fixed_list[i])
            self.move_progress_bar_forward()
            # Update the new constants and reoptimize
            self.get_constants(min_dry_diameter, min_dry_diameter_asymptote, max_dry_diameter_asymptote)
            self.move_progress_bar_forward()
            self.fit_sigmoid_line_to_run()
            self.move_progress_bar_forward()
            # Remake the graph
            figure = self.complete_dry_diameter_graph
            self.create_complete_sigmoid_graph(figure)

            # Update the old processed_data
            self.min_dp_list[self.current_scan] = self.min_dp
            self.min_dp_asym_list[self.current_scan] = self.min_dp_asym
            self.max_dp_asym_list[self.current_scan] = self.max_dp_asym
            self.b_list[self.current_scan] = self.b
            self.d_list[self.current_scan] = self.d
            self.c_list[self.current_scan] = self.c
            self.dp50_list[self.current_scan] = (self.d, self.ss_list[self.current_scan])
            self.updateGUI()
            self.move_progress_bar_forward(complete=True)
        except:
            # Disable the peak
            self.usable_for_kappa_cal_list[self.current_scan] = False
            self.min_dp_list[self.current_scan] = 0
            self.min_dp_asym_list[self.current_scan] = 0
            self.max_dp_asym_list[self.current_scan] = 0
            self.b_list[self.current_scan] = 0
            self.d_list[self.current_scan] = 0
            self.c_list[self.current_scan] = 0
            self.dp50_list[self.current_scan] = (0, 0)
            self.updateGUI()
            self.move_progress_bar_forward(complete=True)

    #----------------------Misc------------------------------------------

    def run(self):
        """start to process the processed_data"""
        self.current_phase = 0
        self.parse_and_match_smps_ccnc_data()

    def reset(self):
        i = 1
        gc.collect()
        for aFigure in self.temperature_graph_list:
            aFigure.clf()
            plt.close(aFigure)
        self.temperature_graph_list = []
        for aFigure in self.adjusted_graph_list:
            plt.close(aFigure)
        self.adjusted_graph_list = []
        for aFigure in self.dry_diameter_graph_list:
            plt.close(aFigure)
        self.dry_diameter_graph_list = []
        if self.min_compare_grpah:
            self.min_compare_grpah.clf()
            plt.close(self.min_compare_grpah)
            self.min_compare_grpah = None
        if self.kappa_graph:
            self.kappa_graph.clf()
            plt.close(self.kappa_graph)
            self.kappa_graph = None
        self.current_point = None
        self.finish_sigmoid_fit_phase = False
        self.cancelling = False
        self.concentration = 0.3
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
        self.kappa_exclude_list = []
        self.kappa_points = []
        gc.collect()

    def move_progress_bar_forward(self, message=None, max_value=None, complete=False, value=1):
        self.view.move_progress_bar_forward(message, max_value, complete, value)
        if self.cancelling is True:
            self.cancelling = False
            raise InterruptError()

    def cancel_current_progress(self):
        self.cancelling = True


def main():
    # controller
    controller = Controller(False)
    #view
    view = View(controller)
    controller.view = view
    view.show_ui()

if __name__ == '__main__':
    main()

