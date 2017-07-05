<<<<<<< HEAD
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
        # The view associated with the controller
        self.view = None
        # The files which contains the data.
        self.files = None
        # The list of paths to the txt files which contain the smps data
        self.smps_txt_file = []
        # The list of paths to the cccn files which contain the ccnc data
        self.ccnc_csv_file = []
        # The peak we are analysing
        self.currPeak = 0
        # The total number of peaks
        self.number_of_peak = 0
        # The single smps file which contains the data
        self.smpsFile = None
        # The single ccnc file which contains the data, created by combining all the ccn files
        self.ccncFile = None
        # The list of start time entries
        self.start_time_list = None
        # the list of end time entries
        self.end_time_list = None
        # The date of the experiment
        self.date = None
        # the time frame of aParam aParam measurement
        self.scan_time = 0
        # The content of the csv file
        self.ccnc_data = None
        # The content of the txt file
        self.smps_data = None
        # The list of all diameter of measurement
        self.diameterList = None
        # The dnlog data
        self.dnlog_list = None
        # The super saturation
        self.ssList = []
        # The SMPS and CCNC data
        self.data = None
        # The raw SMPS and CCNC data
        self.rawData = None
        # The number of process finished
        self.completedStep = 0
        # The current point on the minimum graph
        self.currentPoint = None
        # If the peak data has been optimized
        self.optimized = False
        # If the processing is cancelled
        self.cancel = False

        # Variable needed for peak alignment
        self.flowRate = 0.3

        # variables to store the data
        self.ccnList = None
        self.cnList = None
        self.ccnFixedList = None
        self.cnFixedList = None
        self.gCnList = None
        self.gCcnList = None
        self.ccncSigList = []
        self.ccnNormalizedList = []
        self.diameterMidpointList = []
        self.ccncnSimList = []
        self.temp1 = []
        self.temp2 = []
        self.temp3 = []

        # Variables for peak analysis/optimization/sigmodal fit
        self.b = 0
        self.d = 0
        self.c = 0
        self.minCcn = 4
        self.minDp = 0
        self.maxDp = 0
        self.maxDpAsym = 0
        self.minDpAsym = 0
        self.dp50 = 0
        self.superSaturation = 0
        self.dp50LessCount = 0
        self.dp50MoreCount = 0
        self.dp50Wet = 0
        self.dp50Plus20 = 0
        self.minDpList = []
        self.minDpAsymList = []
        self.maxDpAsymList = []
        self.bList = []
        self.dList = []
        self.cList = []
        self.dp50List = []
        self.dp50WetList = []
        self.dp50Plus20List = []

        # Vars for graphs
        self.adjustedGraph = None
        self.dryDiaGraph = None
        self.minCompareGraph = None
        self.kappaGraph = None
        self.tempGraph = None
        self.adjustedGraphList = []
        self.dryDiaGraphList = []
        self.minPosSMPSList = []
        self.minPosCCNCList = []
        self.tempGraphList = []

        # Vars for kappa
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
        self.kappaExcel = None
        self.shiftList = []
        self.kappaExcludeList = []
        self.kappaPoints = []
        self.klines = None
        self.maxKappa = None
        self.minKappa = None
        self.isFullKappaGraph = False

        # Dicts to store kappa values
        self.kappaCalculatedDict = {}
        self.alphaPineneDict = {}

        # Store the index of the current peak in the data - the first point of the peak
        self.peakPositionInData = []
        # State of a peak, whether usable to calculate kappa or not
        self.usableForKappaCalList = []
        # a clone of self
        self.clone = None


    # -------------General data processing------------------

    def get_flow_rate(self):
        """
        Get the flow rate of the data from user input
        """
        self.flowRate = self.view.InputFlowRate()
        self.flowRate = 1/(self.flowRate*1000/60)

    def get_raw_data(self):
        """
        Get the smps and ccnc data from the input files. 
        """
        self.smps_txt_file = []
        self.ccnc_csv_file = []
        for a_file in self.files:
            if a_file.lower().endswith('.txt'):
                self.smps_txt_file.append(a_file)
            elif a_file.lower().endswith('.csv'):
                self.ccnc_csv_file.append(a_file)
        # check if there are both smps and ccnc files
        if (len(self.smps_txt_file) == 0) or (len(self.ccnc_csv_file) == 0):
            raise FileNotFoundError()
        # if there are, process the file paths, convert to string.
        else:
            self.smps_txt_file = [str(x) for x in self.smps_txt_file]
            self.smps_txt_file = self.smps_txt_file[0]
            self.ccnc_csv_file = [str(x) for x in self.ccnc_csv_file]
    
        # Get the raw contents
        self.ccnc_data = csvProcessing(self.ccnc_csv_file)
        self.date = self.ccnc_data[0]
        self.smps_data = txtProcessing(self.smps_txt_file)
        self.start_time_list = self.smps_data[0]
        # Get the running time (in seconds) of each run
        for i in range(len(self.smps_data)):
            if ''.join(self.smps_data[i][0].split()).lower() == "scanuptime(s)":
                scanUpTime = self.smps_data[i][1]
                scanDownTime = self.smps_data[i+1][1]
                break
        self.scan_time = int(scanUpTime) + int(scanDownTime)
        self.end_time_list = []
        for i in range(len(self.start_time_list)):
            self.end_time_list.append(datetime.strptime(self.start_time_list[i], "%H:%M:%S") + timedelta(
                seconds=self.scan_time))
        #remove the time frames that doesn't have data in the CCNC
        while (True):
            smpsEndTime = self.end_time_list[-1]
            csvEndTime = datetime.strptime(self.ccnc_data[1][len(self.ccnc_data[1]) - 1][0], "%H:%M:%S")
            if smpsEndTime > csvEndTime:
                self.start_time_list = self.start_time_list[:-1]
                self.end_time_list = self.end_time_list[:-1]
            else:
                break
        self.end_time_list = [datetime.strftime(x,"%H:%M:%S") for x in self.end_time_list]
        # Get the total number of run
        self.number_of_peak = len(self.start_time_list)

    def get_dnlog_data(self):
        """
        Get the dnlog data from the txt file
        """
        for k in range(1, len(self.smps_data)):
            if re.search('[aParam-zA-Z]', self.smps_data[k][0]):
                dnlogPos = k
                break
        startTime = ['dp'] + self.start_time_list
        endTime = ['dp'] + self.end_time_list
        dnlog_list = [startTime, endTime]
        for k in range(1, dnlogPos):
            dnlog_list.append(self.smps_data[k][:len(startTime)])
        self.dnlog_list = dnlog_list

    def get_smps_and_ccnc_data(self):
        """
        Get the SMPS data from raw SMPS and CCN from raw CCN
        """
        # Get the data from the SMPS file
        ccnc_data = self.ccnc_data[1]
        smps_data = self.smps_data
        startTime = self.start_time_list
        endTime = self.end_time_list
        width = self.number_of_peak

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

        # Get the CCNC data
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
        self.peakPositionInData = []
        while (True):
            sizeSum = 0
            countSum = 0
            previousTimeStamp = datetime.strptime(startTime[i], "%H:%M:%S").time()
            maxSS = 0
            # Get all data within the time frame
            self.peakPositionInData.append(len(smpsCcnList))
            for t in range(self.scan_time):
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
            k += self.scan_time

            # Update the lists
            self.ssList.append(maxSS)

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
                self.start_time_list = startTime
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
        self.rawData = pandas.DataFrame(smpsCcnList, columns=title)

    def match_smps_ccnc_data(self):
        """
        Calculate the stretch factor of the ccnc measurement
        """
        self.usableForKappaCalList = []
        self.minPosCCNCList = []
        self.minPosSMPSList = []
        additionalDataCount = int(self.scan_time / 2)
        startTime = 0
        endTime = self.scan_time
        currPeak = 0
        newData = [self.rawData.columns.values.tolist()]
        minDist = self.scan_time / 10
        shiftFactor = 0
        self.shiftList = []
        self.number_of_peak = 0
        while True:
            # count currPeak smps
            aPeak = numpy.asarray(self.rawData[startTime:endTime]["SMPS Count"])
            minPosSMPS = getMinIndex(aPeak)
            # assuming that count currPeak of smps is always correct, the first currPeak is in the right position
            if minPosSMPS == -1:
                self.minPosSMPSList.append(None)
                self.minPosCCNCList.append(None)
                self.usableForKappaCalList.append(False)
                minPosSMPS = 0
                minPosCCNC = 0
            else:
                self.minPosSMPSList.append(minPosSMPS + startTime)
                aPeak = numpy.asarray(self.rawData[startTime + shiftFactor:endTime + additionalDataCount]["CCNC Count"])
                minPosCCNC = getMinIndex(aPeak)
                if minPosCCNC == -1:
                    self.minPosSMPSList[-1] = None
                    self.minPosCCNCList.append(None)
                    self.usableForKappaCalList.append(False)
                    minPosSMPS = 0
                    minPosCCNC = 0
                else:
                    self.usableForKappaCalList.append(True)
                    shiftFactor = shiftFactor + minPosCCNC - minPosSMPS
                    if currPeak == 0:
                        self.minPosCCNCList.append(minPosSMPS)
                    else:
                        self.minPosCCNCList.append(minPosCCNC + minPosCCNC - minPosSMPS + startTime)
            self.shiftList.append(shiftFactor)

            # Add all peaks to data, whether a peak is a good one or not.
            for i in range(self.scan_time):
                scanTime = i + 1
                realTime = self.rawData.iat[startTime + i, 1]
                dp = self.rawData.iat[startTime + i, 2]
                smpsCount = self.rawData.iat[startTime + i, 3]
                ccncCount = self.rawData.iat[startTime + i + shiftFactor, 4]
                aveSize = self.rawData.iat[startTime + i + shiftFactor, 5]
                t1 = self.rawData.iat[startTime + i + shiftFactor, 6]
                t2 = self.rawData.iat[startTime + i + shiftFactor, 7]
                t3 = self.rawData.iat[startTime + i + shiftFactor, 8]
                newData.append([scanTime, realTime, dp, smpsCount, ccncCount, aveSize,t1,t2,t3])

            currPeak += 1
            self.number_of_peak += 1
            if endTime + self.scan_time + shiftFactor >= len(self.rawData) or currPeak >= len(self.peakPositionInData):
                break

            startTime = self.peakPositionInData[currPeak]
            endTime = startTime + self.scan_time

        headers = newData.pop(0)
        self.data = pandas.DataFrame(newData, columns=headers)
        self.data = dateConvert(self.data)

    def peak_align_and_graph(self):
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
        figure.canvas.mpl_connect('pick_event', self.onPick)

        # Prepare the data
        tempSMPSPeakList = []
        tempCCNPeakCList = []
        # Plot only valid points
        for i in range(len(self.minPosCCNCList)):
            if self.minPosCCNCList[i] and self.minPosSMPSList[i]:
                tempSMPSPeakList.append(self.minPosSMPSList[i])
                tempCCNPeakCList.append(self.minPosCCNCList[i])
            else:
                # Make up for the null values
                if len(tempSMPSPeakList) >0:
                    tempSMPSPeakList.append(tempSMPSPeakList[-1] + self.scan_time)
                    tempCCNPeakCList.append(tempCCNPeakCList[-1] + self.scan_time)
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
        self.currentPoint, = plt.plot(x[0],y[0],'o', color = "#81C784", ms = 12, mew = 0)
        plt.xlabel("SMPS minumum point")
        plt.ylabel("CCNC minimum point")
        handles, labels = plt.gca().get_legend_handles_labels()
        legend = plt.legend(handles, labels,loc="upper right", bbox_to_anchor=(1,0.7))
        legend.get_frame().set_facecolor('#9E9E9E')

        if not self.minCompareGraph:
            self.minCompareGraph = plt.gcf()
        else:
            plt.close(self.minCompareGraph)
            self.minCompareGraph = plt.gcf()

    def onPick(self, event):
        """
        When aParam dot on the graph is selected
        """
        peak = event.ind[0]
        if peak != self.currPeak:
            self.switchToPeak(peak)

    # -------------Single peak processing ------------------

    def preparePeakData(self):
        """
        Convert data to separate lists for easier access. Also normalize and calculate necessary lists. Prepare data
        for later usage
        """
        try:
            ccnCNList = []
            dNdLogDpList = []
            startPoint = self.currPeak * self.scan_time
            endPoint = (self.currPeak + 1) * self.scan_time
            self.diameterList = list(self.data[startPoint:endPoint]['dp'])
            ccnList = list(self.data[startPoint:endPoint]['CCNC Count'])
            cnList = list(self.data[startPoint:endPoint]['SMPS Count'])
            self.temp1 = list(self.data[startPoint:endPoint]['T1'])
            self.temp2 = list(self.data[startPoint:endPoint]['T2'])
            self.temp3 = list(self.data[startPoint:endPoint]['T3'])

            # modify this to get aParam decently good number
            cnList = [x * 0.2 for x in cnList]
            self.ccnList = ccnList
            self.cnList = cnList
            self.diameterMidpointList = []
            self.ccncnList = ccnCNList
            for i in range(len(ccnList)):
                self.ccncnList.append(ccnList[i] / cnList[i])
            self.dropSizeList = list(self.data[startPoint:endPoint]['Ave Size'])
            for i in range(2, len(self.dnlog_list)):
                self.diameterMidpointList.append(self.dnlog_list[i][0])
                dNdLogDpList.append(self.dnlog_list[i][self.currPeak + 1])
            self.ccnNormalizedList = normalizeList(dNdLogDpList)
            self.usableForKappaCalList[self.currPeak] = True
        except:
            raise OptimizationError()

    def initCorrectCharges(self):
        """
        Initiate the correct charge procedure
        """
        # If the peak doesn't have good alignment, do nothing
        self.cnFixedList = self.cnList[:]
        self.ccnFixedList = self.ccnList[:]
        self.gCcnList = self.ccnFixedList[:]
        self.gCnList = self.cnFixedList[:]

    def correctCharges(self):
        """
        Perform charge correction
        """
        try:
            # If the peak doesn't have good alignment, do nothing
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
            # frac0List = fractionCalculation(self.diameterList,0,coeficientList[0])
            frac1List = fractionCalculation(self.diameterList, 1, coeficientList[1])
            frac2List = fractionCalculation(self.diameterList, 2, coeficientList[2])
            frac3List = fractionCalculation(self.diameterList, 3)
            chargeList = []

            for i in self.diameterList:
                QtGui.qApp.processEvents()
                aDList = [0]
                for k in range(1, 4):
                    c = calCC(i * 10 ** -9, lambdaAir)
                    dp = 10 ** 9 * findDp(i * 10 ** -9 / c, lambdaAir, k)
                    aDList.append(dp)
                chargeList.append(aDList)

            # second part of correct charges
            self.cnFixedList = self.cnList[:]
            self.ccnFixedList = self.ccnList[:]
            maxUpperBinBound = (self.diameterList[-1] + self.diameterList[-2]) / 2
            lenDpList = len(self.diameterList)

            for i in range(lenDpList):
                QtGui.qApp.processEvents()
                n = lenDpList - i - 1
                moveDoubletCounts = frac2List[n] / (frac1List[n] + frac2List[n] + frac3List[n]) * self.cnList[n]
                moveTripletCounts = frac3List[n] / (frac1List[n] + frac2List[n] + frac3List[n]) * self.cnList[n]
                self.cnFixedList[n] = self.cnFixedList[n] - moveDoubletCounts - moveTripletCounts
                self.ccnFixedList[n] = self.ccnFixedList[n] - moveDoubletCounts - moveTripletCounts
                if chargeList[n][2] <= maxUpperBinBound:
                    j = lenDpList - 2
                    while (True):
                        upperBinBound = (self.diameterList[j] + self.diameterList[j + 1]) / 2
                        lowerBinBound = (self.diameterList[j] + self.diameterList[j - 1]) / 2
                        if upperBinBound > chargeList[n][2] >= lowerBinBound:
                            self.cnFixedList[j] = self.cnFixedList[j] + moveDoubletCounts
                            if chargeList[n][2] < asymp:
                                if self.gCcnList[j] > epsilon:
                                    self.ccnFixedList[j] = self.ccnFixedList[j] + moveDoubletCounts * self.gCcnList[j] / \
                                                                                  self.gCnList[j]
                            else:
                                self.ccnFixedList[j] = self.ccnFixedList[j] + moveDoubletCounts
                            break
                        j -= 1
                if chargeList[n][3] < maxUpperBinBound:
                    j = lenDpList - 2
                    while (True):
                        upperBinBound = (self.diameterList[j] + self.diameterList[j + 1]) / 2
                        lowerBinBound = (self.diameterList[j] + self.diameterList[j - 1]) / 2
                        if upperBinBound > chargeList[n][3] >= lowerBinBound:
                            self.cnFixedList[j] = self.cnFixedList[j] + moveTripletCounts
                            if chargeList[n][3] < asymp:
                                self.ccnFixedList[j] = self.ccnFixedList[j] + moveTripletCounts * self.ccnList[j] / \
                                                                              self.cnList[j]
                            else:
                                self.ccnFixedList[j] = self.ccnFixedList[j] + moveTripletCounts
                            break
                        j -= 1

            for i in range(len(self.ccnFixedList)):
                QtGui.qApp.processEvents()
                if self.ccnFixedList[i] / self.cnFixedList[i] < -0.01:
                    self.ccnFixedList[i] = 0

            self.gCcnList = self.ccnFixedList[:]
            self.gCnList = self.cnFixedList[:]
            self.usableForKappaCalList[self.currPeak] = True
        except:
            raise OptimizationError()

    def getConstants(self, minDp = None, minDPAsym = None, maxDpAsym = None):
        """
        Acquire the necessary constants from the data to perform sigmodal fit
        """
        try:
            if minDp and minDPAsym and maxDpAsym:
                self.minDp = minDp
                self.minDpAsym = minDPAsym
                self.maxDpAsym = maxDpAsym
            else:
                asymList = getAsym(self.diameterList, self.ccncSigList)
                mpIndex = 0

                # Find the index of midPoint
                checkLength = self.scan_time / 20
                for i in range(checkLength,len(self.ccncSigList)-checkLength):
                    isMid = False
                    if self.ccncSigList[i] > 0.5:
                        isMid = True
                        # check if the previous 5 numbers are smaller and the next 5 numbers are bigger
                        for j in range(1,checkLength):
                            if self.ccncSigList[i+j] < self.ccncSigList[i]:
                                isMid = False
                                break
                    if isMid:
                        mpIndex = i
                        break

                minDpAsymPos = 0
                currMax = 0
                # Get minDP
                for i in range(mpIndex,1, -1):
                    if self.ccncSigList[i] < 0.1:
                        self.minDp = self.diameterList[i]
                        break

                # Get minDpAsym
                for i in range(mpIndex, mpIndex + self.scan_time / 10):
                    if i < len(self.diameterList):
                        if self.ccncSigList[i] > currMax:
                            minDpAsymPos = i
                            currMax = self.ccncSigList[i]
                            self.minDpAsym = self.diameterList[i]

                maxCCN = max(self.ccncSigList)
                if maxCCN <= 1.2:
                    stableThreshold = 0.075
                else:
                    stableThreshold = 0.25
                # determine maxDplen(self.ccncSigList)
                for i in range(minDpAsymPos + 5, len(self.ccncSigList)):
                    if self.ccncSigList[i] > 1.3:
                        self.maxDpAsym = self.diameterList[i]
                        break
                    elif abs(self.ccncSigList[i] - self.ccncSigList[i - 1]) > 2 * stableThreshold:
                        self.maxDpAsym = self.diameterList[i - 1]
                        break

            self.maxDp = self.maxDpAsym
            # Get the data
            asymsList = []
            for i in range(len(self.diameterList)):
                if self.minDpAsym < self.diameterList[i] < self.maxDpAsym:
                    asymsList.append(self.ccncSigList[i])
                else:
                    asymsList.append(0)
            # Calcualte constants
            self.b = getAveNoneZero(asymsList)
            self.ccncnSimList.append(0)
            for i in range(1, len(self.diameterList)):
                if self.minDp < self.diameterList[i] < self.maxDp:
                    n = self.b / (1 + (self.diameterList[i] / self.d) ** self.c)
                    self.ccncnSimList.append(n)
                else:
                    self.ccncnSimList.append(self.ccncnSimList[i - 1])
            self.usableForKappaCalList[self.currPeak] = True
        except:
            raise OptimizationError()

    def optimizePeak(self):
        """
        Perform optimization, which is basically sigmodal fit
        """
        try:
            xList = []
            yList = []
            for i in range(len(self.diameterList)):
                if self.minDp < self.diameterList[i] < self.maxDp:
                    xList.append(self.diameterList[i])
                    yList.append(self.ccncSigList[i])
            initGuess = [30, -10]
            for i in range(len(yList)):
                if math.isnan(yList[i]) or math.isinf(yList[i]):
                    yList[i] = 0
            xList = numpy.asarray(xList)
            yList = numpy.asarray(yList)
            result = opt.curve_fit(f, xList, yList, bounds=([self.minDp, -200], [self.maxDp, -1]), method="trf")
            self.d = result[0][0]
            self.c = result[0][1]
            self.ccncnSimList = [0]
            for i in range(1, len(self.diameterList)):
                if self.diameterList[i] > self.d:
                    self.dp50Wet = self.dropSizeList[i - 1]
                    break
            for i in range(1, len(self.diameterList)):
                if self.diameterList[i] > (self.d + 20):
                    self.dp50Plus20 = self.dropSizeList[i - 1]
                    break

            for i in range(1, len(self.diameterList)):
                if self.minDp < self.diameterList[i] < self.maxDp:
                    n = self.b / (1 + (self.diameterList[i] / self.d) ** self.c)
                    self.ccncnSimList.append(n)
                else:
                    self.ccncnSimList.append(self.ccncnSimList[i - 1])
            self.usableForKappaCalList[self.currPeak] = True
        except:
            raise OptimizationError()

    # -------------Produce graphs ------------------

    def makeAdjustedGraph(self, newFigure = None):
        """
        Make the comparable CCN and CN graph after adjustment of minimum peak
        """
        try:
            if len(self.cnList) != len(self.ccnList) or len(self.cnList) == 0 or len(self.ccnList) == 0:
                return
            if newFigure is None:
                figure = plt.figure(facecolor=settings.graphBackgroundColor)
            else:
                figure = newFigure
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
            x = range(self.scan_time)
            plt.plot(x, self.cnList, linewidth = 4, color = '#EF5350', label = "CN")
            plt.plot(x, self.ccnList, linewidth = 4, color = '#2196F3',label = "CCN")
            handles, labels = plt.gca().get_legend_handles_labels()
            legend = plt.legend(handles, labels, loc="upper left", bbox_to_anchor=(0, 0.9))
            legend.get_frame().set_facecolor('#9E9E9E')

            plt.xlabel("Scan time(s)")
            plt.ylabel("Concentration(cm3)")
        except:
            figure = plt.figure(facecolor=settings.graphBackgroundColor)
        finally:
            self.adjustedGraphList.append(plt.gcf())

    def makeCCNGraph(self, newFigure = None):
        """
        Make graph for the CCN data. Used to cross check for peak alignment.
        """
        try:
            if newFigure is None:
                figure = plt.figure(facecolor=settings.graphBackgroundColor)
            else:
                figure = newFigure
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
            plt.plot(self.diameterList, self.ccncnList, 'o', color="#2196F3", mew=0.5,mec = "#0D47A1",  ms = 9, label = "CCN/CN")
            yLim = min(2, max(self.ccncnList)) + 0.2
            plt.gca().axes.set_ylim([-0.1, yLim])
            handles, labels = plt.gca().get_legend_handles_labels()
            legend = plt.legend(handles, labels,loc="upper right", bbox_to_anchor=(1,0.9))
            legend.get_frame().set_facecolor('#9E9E9E')
            plt.xlabel("Diameter (nm)")
            plt.ylabel("CCN/CN")
        except:
            figure = plt.figure(facecolor=settings.graphBackgroundColor)
        finally:
            self.dryDiaGraphList.append(plt.gcf())

    def makeTempGraph(self, newFigure = None):
        """
        Make the temperature graphs
        """
        try:
            if newFigure is None:
                figure = plt.figure(facecolor=settings.graphBackgroundColor)
            else:
                figure = newFigure
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

            x = range(self.scan_time)
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
            self.tempGraphList.append(plt.gcf())

    def makeFullDryDiameterGraph(self, newFigure = None):
        """
        Make complete graph of the dry diameter after optimization and sigmodal fit
        """
        try:
            if not self.minPosCCNCList[self.currPeak] or not self.minPosSMPSList[self.currPeak]:
                self.makeCCNGraph(newFigure)
                return
            if newFigure is None:
                figure = plt.figure(facecolor=settings.graphBackgroundColor)
            else:
                figure = newFigure
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
            yLim = min(2, max(self.ccncSigList)) + 0.2
            plt.gca().axes.set_ylim([-0.1, yLim])

            plt.plot(self.diameterMidpointList, self.ccnNormalizedList, linewidth=4, color='#43A047', label="dN/dLogDP")
            # Graph the sig fit only if the peak is usable
            if self.usableForKappaCalList[self.currPeak]:
                plt.plot(self.diameterList, self.ccncnSimList, linewidth=5, color='#EF5350', label="Sigmodal Fit")
            plt.plot(self.diameterList, self.ccncnList, 'o', color="#2196F3", mew=0.5,mec = "#1976D2",  ms = 9, label = "CCN/CN")
            plt.plot(self.diameterList, self.ccncSigList, 'o', color="#1565C0", mew=0.5,mec = "#0D47A1",  ms = 9, label = "CCN/CN (Corrected)")
            plt.xlabel("Dry diameter(nm)")
            plt.ylabel("CCN/CN ratio and Normalized dN/dLogDP")
            handles, labels = plt.gca().get_legend_handles_labels()
            legend = plt.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.7, 1.1))
            legend.get_frame().set_facecolor('#9E9E9E')
            self.dryDiaGraphList.append(plt.gcf())
        except:
            figure = plt.figure(facecolor=settings.graphBackgroundColor)
            self.dryDiaGraphList.append(plt.gcf())

    def makeInitialGraphs(self):
        """
        Make initial graph without optimizations
        """
        plt.ioff()
        self.makeAdjustedGraph()
        self.makeCCNGraph()
        self.makeTempGraph()

    def makeOptimizedGraphs(self):
        """
        Make complete graphs after optimization
        """
        plt.ioff()
        self.makeFullDryDiameterGraph()

    def makeKappaGraph(self):
        """
        Produce the kappa graph, may be in full or only around the points
        """
        # Read in the kappa lines csv files.
        if self.klines is None:
            self.klines = pandas.read_csv("klines.csv", header=1)
        # Acquire the kappa points
        kappaList = []
        stdKappaList = []
        self.kappaPoints = []

        for aKey in self.alphaPineneDict.keys():
            kappaList.append(self.alphaPineneDict[aKey][2])
            stdKappaList.append(self.alphaPineneDict[aKey][3])
            if not self.isFullKappaGraph:
                self.kappaPoints.append((self.alphaPineneDict[aKey][0], aKey))
            else:
                for i in self.alphaPineneDict[aKey][-1]:
                    self.kappaPoints.append((i, aKey))

        # Determine maximum and minimum kappas
        tempKList = []
        for i in range(len(kappaList)):
            tempKList.append(kappaList[i] + stdKappaList[i])
        self.maxKappa = max(tempKList)
        tempKList = []
        for i in range(len(kappaList)):
            tempKList.append(kappaList[i] - stdKappaList[i])
        self.minKappa = min(tempKList)
        header = self.klines.columns
        diaList = self.klines[header[1]]

        fullKXList = []
        fullKYList = []
        kpXList = []
        kpYList = []
        excludedXList = []
        excludedYList = []

        for i in range(len(self.kappaPoints)):
            fullKXList.append(self.kappaPoints[i][0])
            fullKYList.append(self.kappaPoints[i][1])
        # Process which point is included and which point is excluded
        # Include all excluded points to exclude list
        for i in self.kappaExcludeList:
            excludedXList.append(self.kappaPoints[i][0])
            excludedYList.append(self.kappaPoints[i][1])
        # Exclude all excluded points from kpx/y list
        for i in range(len(self.kappaPoints)):
            if i not in self.kappaExcludeList:
                kpXList.append(self.kappaPoints[i][0])
                kpYList.append(self.kappaPoints[i][1])
        # convert list to array
        excludedXList = numpy.asarray(excludedXList)
        excludedYList = numpy.asarray(excludedYList)
        kpXList = numpy.asarray(kpXList)
        kpYList = numpy.asarray(kpYList)

        # Get all kappa points - without average
        # Prepare the figure
        if self.kappaGraph:
            figure = plt.figure(self.kappaGraph.number)
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
        figure.canvas.mpl_connect('pick_event', self.onKappaPick)

        # Graph the klines
        # Full kappa lines
        #if self.isFullKappaGraph:
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
            if self.maxKappa > kappa:
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
            if self.minKappa > kappa:
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
        self.kappaGraph = plt.gcf()
        self.view.updateKappaGraph()

    def onKappaPick(self, event):
        """
        When a kappa point is clicked
        """
        kappaPoint = event.ind[0]
        excluded = False
        # if already in excluded list, include the points
        for i in range(len(self.kappaExcludeList)):
            if self.kappaExcludeList[i] == kappaPoint:
                self.kappaExcludeList = self .kappaExcludeList[:i] + self.kappaExcludeList[i+1:]
                excluded = True
                break
        # else, exclude the point
        if not excluded:
            self.kappaExcludeList.append(kappaPoint)
        # Remake the kappa graph
        self.makeKappaGraph()

    # -------------Single peak procedure ------------------

    def singlePeakProcessingProcedure(self):
        """
        The steps required to process a single peak
        """
        try:
            self.preparePeakData()
            self.makeProgress()
            self.makeInitialGraphs()
            self.makeProgress()
        except OptimizationError:
            # Make the peak invalid
            self.minPosCCNCList[self.currPeak] = None
            self.minPosSMPSList[self.currPeak] = None
            self.usableForKappaCalList[self.currPeak] = False

    def singlePeakOptimizationProcedure(self):
        """
        The steps required to optimize a single peak
        """
        try:
            # If a peak is invalid, then do nothing
            if not self.minPosCCNCList[self.currPeak] or not self.minPosSMPSList[self.currPeak]:
                self.preparePeakData()
                raise OptimizationError()
            else:
                self.ccncSigList = []
                self.preparePeakData()
                self.makeProgress()
                removeSmallCcn(self.ccnList, self.minCcn)
                self.makeProgress()
                self.initCorrectCharges()
                for i in range(5):
                    self.correctCharges()
                    self.makeProgress()
                for i in range(len(self.ccnFixedList)):
                    self.ccncSigList.append(self.ccnFixedList[i] / self.cnFixedList[i])
                self.makeProgress()
                self.getConstants()
                self.makeProgress()
                self.optimizePeak()
                self.makeProgress()
                self.makeOptimizedGraphs()
                self.makeProgress()
        except OptimizationError:
            # Disable the peak
            self.usableForKappaCalList[self.currPeak] = False
            # Add empty data to the list
            self.bList.append(0)
            self.dList.append(0)
            self.cList.append(0)
            self.dp50List.append((0, 0))
            self.dp50WetList.append(0)
            self.dp50Plus20List.append(0)
            self.minDpList.append(0)
            self.minDpAsymList.append(0)
            self.maxDpAsymList.append(0)
            self.makeFullDryDiameterGraph()
        else:
            # Store data
            self.minDpList.append(self.minDp)
            self.minDpAsymList.append(self.minDpAsym)
            self.maxDpAsymList.append(self.maxDpAsym)
            self.bList.append(self.b)
            self.dList.append(self.d)
            self.cList.append(self.c)
            self.dp50WetList.append(self.dp50Wet)
            self.dp50Plus20List.append(self.dp50Plus20)
            self.dp50List.append((self.d, self.ssList[self.currPeak]))

    # --------------General processing procedure-------------

    def preprocess_data(self):
        """
        Pre-process the data. This includes extract the dnlog data, smps and ccnc data
        from the raw data file. 
        """
        try:
            # reset all constants and the GUI
            self.reset()
            self.view.reset()
            #initiate the progress bar
            self.makeProgress(maxValue=5)
            self.get_flow_rate()
            self.makeProgress("Processing raw data from files...")
            self.get_raw_data()
            self.makeProgress("Acquiring dnlog data...")
            self.get_dnlog_data()
            self.makeProgress("Acquiring SMPS and CCNC data...")
            self.get_smps_and_ccnc_data()
            self.makeProgress("Transforming the CCNC data to match SMPS data....")
            self.match_smps_ccnc_data()
            self.completedStep = 1
            self.makeProgress(complete=True)
        except FileNotFoundError:
            self.view.showError("Can't find SMPS or CCNC data in files: " + str(self.files))
            raise DataPreparationError()
        except FileProcessingError:
            self.view.showError("Can't process the SMPS or CCNC files!")
            raise DataPreparationError()
        except dnlogDataError:
            self.view.showError("Can't process dnlog data from the SMPS file!")
            raise DataPreparationError()
        except DataError:
            self.view.showError("Can't process data from SMPS or CCNC raw data!")
            raise DataPreparationError()
        except DataMatchingError:
            self.view.showError("Can't match SMPS and CCNC data!")
            raise DataPreparationError()
        except InterruptError:
            self.reset()
            self.view.reset()
            self.view.showError("The preparation process is cancelled")
            raise DataPreparationError()

    def initialProcessProcedure(self):
        """
        The processes required to process all the peak data
        """
        try:
            if self.completedStep < 1:
                self.view.showError("Data preparation process is not completed.")
                raise InterruptError
            self.makeProgress("Preparing data for peak processing...", maxValue=self.number_of_peak * 3)
            for i in range(0, self.number_of_peak):
                self.currPeak = i
                self.makeProgress("Processing peak " + str(i + 1), value=0)
                self.singlePeakProcessingProcedure()
            self.currPeak = -1
            self.makeProgress(complete=True)
            self.view.updateGeneralInfo()
            self.peak_align_and_graph()
            self.switchToPeak(0)
            self.completedStep = 2
        except InterruptError:
            self.view.showError("The data initialization process is cancelled")

    def optimizationProcedure(self):
        """
        The processes required to optimize/sigmodal fit all peaks
        """
        try:
            if self.completedStep < 2:
                self.view.showError("Data initialization process is not completed.")
                return

            #clear temp graphs
            for aFigure in self.tempGraphList:
                aFigure.clf()
                plt.close(aFigure)
            self.tempGraphList = []
            # clear dry diameter graphs
            for aFigure in self.dryDiaGraphList:
                aFigure.clf()
                plt.close(aFigure)
            self.dryDiaGraphList = []
            self.makeProgress("Preparing data for peak processing...", maxValue=self.number_of_peak * 12)
            self.optimized = True
            for i in range(0, self.number_of_peak):
                self.currPeak = i
                self.makeProgress("Processing peak " + str(i + 1), value=0)
                self.singlePeakOptimizationProcedure()
            # reset currPeak
            self.currPeak = -1
            self.completedStep = 3
            self.makeProgress(complete=True)
            self.switchToPeak(0)
        except InterruptError:
            self.view.showError("The optimization process is cancelled")
            self.optimized = False
            for aFigure in self.tempGraphList:
                aFigure.clf()
                plt.close(aFigure)
            self.tempGraphList = []
            # clear dry diameter graphs
            for aFigure in self.dryDiaGraphList:
                aFigure.clf()
                plt.close(aFigure)
            self.dryDiaGraphList = []
            self.initialProcessProcedure()



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
        Calculate the kappa values - producing both raw data kappa and graph data kappa
        """
        for i in range(len(self.dp50List)):
            self.usableForKappaCalList.append(True)
        self.makeProgress("Calculating Kappa...", maxValue=len(self.dp50List) + 3)
        self.makeProgress("Reading in lookup table")
        if self.kappaExcel is None:
            self.kappaExcel = pandas.read_csv("kCal.csv", header=None)
        lookup = self.kappaExcel
        self.makeProgress("Calculating basic consts...")
        self.aParam = 0.00000869251 * self.sigma / self.temp
        self.asc = (exp(sqrt(4 * self.aParam ** 3 / (27 * self.iKappa * (self.dd * 0.000000001) ** 3))) - 1) * 100
        # Calculate each kappa
        firstAKappa = 0
        scCalcs = False
        for i in range(len(self.dp50List)):
            if not self.usableForKappaCalList[i]:
                continue
            ss = float(self.dp50List[i][1])
            dp50 = float(self.dp50List[i][0])
            rowIndex = int(math.floor(dp50 - 9))
            matchRow = list(lookup.iloc[rowIndex][1:])
            valueRow = list(lookup.iloc[0][1:])
            a = getCorrectNum(matchRow, ss)
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
                self.makeProgress("Calculating sc...")
                # Calculate the sc Calculation
                firstAKappa = self.appKappa
                # Calcualte the first row of scCalcs
                for i in range(len(self.dp50List)):
                    if self.dp50List[i][0] != 0:
                        dList1 = [float(self.dp50List[i][0]) * 0.000000001]
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
            self.makeProgress("Processing peak" + str(i + 1))
            self.trueSC = (max(sList1[i:]) - 1) * 100
            self.anaKappa = (4 * self.aParam ** 3) / (27 * (dp50 * 0.000000001) ** 3 * log(ss / 100 + 1) ** 2)
            kDevi = (self.appKappa - self.anaKappa) / self.appKappa * 100
            if ss in self.kappaCalculatedDict.keys():
                self.kappaCalculatedDict[ss].append([dp50, self.appKappa, self.anaKappa, kDevi, self.trueSC])
            else:
                self.kappaCalculatedDict[ss] = ([[dp50, self.appKappa, self.anaKappa, kDevi, self.trueSC]])

        for aKey in self.kappaCalculatedDict.keys():
            aSSList = self.kappaCalculatedDict[aKey]
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
            self.alphaPineneDict[aKey] = (meanDp, stdDp, meanApp, stdApp, meanAna, stdAna, meanDev, devMean,dp50List)
        self.makeProgress(complete=True)

    #--------------- UI interaction - Data interaction functions--------------------

    def updateGUI(self):
        """
        Update the GUI with the newest data of the current peak
        """
        self.adjustedGraph = self.adjustedGraphList[self.currPeak]
        self.dryDiaGraph = self.dryDiaGraphList[self.currPeak]
        self.view.update_dp_dnlog_figures(self.adjustedGraph, self.dryDiaGraph)
        if self.optimized:
            self.b = self.bList[self.currPeak]
            self.d = self.dList[self.currPeak]
            self.c = self.cList[self.currPeak]
            self.minDpAsym = self.minDpAsymList[self.currPeak]
            self.minDp = self.minDpList[self.currPeak]
            self.maxDpAsym = self.maxDpAsymList[self.currPeak]
            self.superSaturation = self.dp50List[self.currPeak][1]
            self.dp50 = self.dp50List[self.currPeak][0]
            self.dp50Wet = self.dp50WetList[self.currPeak]
            self.dp50Plus20 = self.dp50Plus20List[self.currPeak]
            self.dp50LessCount = 0
            self.dp50MoreCount = 0
            for i in self.diameterList:
                if i >= self.dp50:
                    self.dp50MoreCount += 1
                else:
                    self.dp50LessCount += 1
            self.currentPoint.set_xdata(numpy.asarray(self.minPosSMPSList)[self.currPeak])
            self.currentPoint.set_ydata(numpy.asarray(self.minPosCCNCList)[self.currPeak])
            self.view.updateTempOrMinFigure(self.minCompareGraph)
            self.view.updateSigFitPeakInfo()

        else:
            self.tempGraph = self.tempGraphList[self.currPeak]
            self.superSaturation = self.ssList[self.currPeak]
            self.view.updateTempOrMinFigure(self.tempGraph)
            self.view.updateBasicPeakInfo()

    def shiftOneSecond(self, forward=True):
        try:
            if self.completedStep <=1:
                return
            # If the peak is invalid, do nothing
            if not self.usableForKappaCalList[self.currPeak] or self.minPosSMPSList[self.currPeak] is None or \
                            self.minPosCCNCList[self.currPeak] is None:
                return
            # Get the shift factor
            shiftFactor = self.shiftList[self.currPeak]
            # Change the shift factor
            if forward == True:
                shiftFactor -= 1
            else:
                shiftFactor += 1
            self.shiftList[self.currPeak] = shiftFactor
            # Get the new data list from raw data
            startTime = self.peakPositionInData[self.currPeak] + shiftFactor
            endTime = startTime + self.scan_time
            newCCN = list(self.rawData.iloc[startTime:endTime, 4])
            newAveSize = list(self.rawData.iloc[startTime:endTime, 5])
            newT1 = list(self.rawData.iloc[startTime:endTime,6])
            newT2 = list(self.rawData.iloc[startTime:endTime,7])
            newT3 = list(self.rawData.iloc[startTime:endTime,8])

            # Get the location of the replacement in the data
            dataStartTime = self.scan_time * self.currPeak
            dataEndTime = dataStartTime + self.scan_time

            # Replace the data in the data
            for i in range(self.scan_time):
                self.data.iat[dataStartTime + i, 4] = newCCN[i]
                self.data.iat[dataStartTime + i, 5] = newAveSize[i]
                self.data.iat[dataStartTime + i, 6] = newT1[i]
                self.data.iat[dataStartTime + i, 7] = newT2[i]
                self.data.iat[dataStartTime + i, 8] = newT3[i]

            # Update the data
            self.preparePeakData()
            # Fix the graph
            figure = self.adjustedGraph
            self.makeAdjustedGraph(newFigure=figure)
            figure = self.dryDiaGraph
            self.makeCCNGraph(figure)
            figure = self.tempGraph
            self.makeTempGraph(figure)
            # Update the graph in view
            self.view.updateDpdnlogFigures(self.adjustedGraph, self.dryDiaGraph)
            self.view.updateTempOrMinFigure(self.tempGraph)
        except:
            self.view.showError("Can't shift data. You should disable this peak!")

    def disablePeak(self):
        # disable the peak
        if self.completedStep == 0:
            return
        if not self.optimized:
            self.minPosSMPSList[self.currPeak] = None
            self.minPosCCNCList[self.currPeak] = None
            self.usableForKappaCalList[self.currPeak] = False
        else:
            self.usableForKappaCalList[self.currPeak] = False
        self.updateGUI()

    def switchToPeak(self, peak):
        """
        Switch to a new peak
        :param peak: the new peak
        """
        if peak != self.currPeak:
            self.currPeak = peak
            self.updateGUI()

    def reOptimization(self, minDp, minDpAsym, maxDpAsym):
        """
        Re-optimize the data after changing constants minDp, minDpAsym, maxDpAsym
        """
        if not self.optimized:
            return
        self.makeProgress("Re-Optimizaing peak" + str(self.currPeak + 1), maxValue=6)
        self.usableForKappaCalList[self.currPeak] = True
        self.ccncSigList = []
        # Get data
        try:
            self.preparePeakData()
            self.makeProgress()
            removeSmallCcn(self.ccnList, self.minCcn)
            self.makeProgress()
            self.initCorrectCharges()
            for i in range(5):
                self.correctCharges()
                self.makeProgress()
            for i in range(len(self.ccnFixedList)):
                self.ccncSigList.append(self.ccnFixedList[i] / self.cnFixedList[i])
            self.makeProgress()
            # Update the new constants and reoptimize
            self.getConstants(minDp, minDpAsym, maxDpAsym)
            self.makeProgress()
            self.optimizePeak()
            self.makeProgress()
            # Remake the graph
            figure = self.dryDiaGraph
            self.makeFullDryDiameterGraph(figure)

            # Update the old data
            self.minDpList[self.currPeak] = self.minDp
            self.minDpAsymList[self.currPeak] = self.minDpAsym
            self.maxDpAsymList[self.currPeak] = self.maxDpAsym
            self.bList[self.currPeak] = self.b
            self.dList[self.currPeak] = self.d
            self.cList[self.currPeak] = self.c
            self.dp50List[self.currPeak] = (self.d, self.ssList[self.currPeak])
            self.updateGUI()
            self.makeProgress(complete=True)
        except:
            # Disable the peak
            self.usableForKappaCalList[self.currPeak] = False
            self.minDpList[self.currPeak] = 0
            self.minDpAsymList[self.currPeak] = 0
            self.maxDpAsymList[self.currPeak] = 0
            self.bList[self.currPeak] = 0
            self.dList[self.currPeak] = 0
            self.cList[self.currPeak] = 0
            self.dp50List[self.currPeak] = (0, 0)
            self.updateGUI()
            self.makeProgress(complete=True)

    #----------------------Misc------------------------------------------

    def run(self):
        """start to process the data"""
        self.completedStep = 0
        try:
            self.preprocess_data()
        except DataPreparationError:
            pass
        else:
            self.initialProcessProcedure()

    def reset(self):
        """
        Reset all values of the program to process a new data set
        Clear out the memory as well
        """
        i = 1
        gc.collect()
        for aFigure in self.tempGraphList:
            aFigure.clf()
            plt.close(aFigure)
        self.tempGraphList = []
        for aFigure in self.adjustedGraphList:
            plt.close(aFigure)
        self.adjustedGraphList = []
        for aFigure in self.dryDiaGraphList:
            plt.close(aFigure)
        self.dryDiaGraphList = []
        if self.minCompareGraph:
            self.minCompareGraph.clf()
            plt.close(self.minCompareGraph)
            self.minCompareGraph = None
        if self.kappaGraph:
            self.kappaGraph.clf()
            plt.close(self.kappaGraph)
            self.kappaGraph = None
        self.currentPoint = None
        self.optimized = False
        self.cancel = False
        self.flowRate = 0.3
        self.ccnList = None
        self.cnList = None
        self.ccnFixedList = None
        self.cnFixedList = None
        self.gCnList = None
        self.gCcnList = None
        self.ccncSigList = []
        self.ccnNormalizedList = []
        self.diameterMidpointList = []
        self.ccncnSimList = []
        self.b = 0
        self.d = 0
        self.c = 0
        self.minCcn = 4
        self.minDp = 0
        self.maxDp = 0
        self.maxDpAsym = 0
        self.minDpAsym = 0
        self.dp50 = 0
        self.superSaturation = 0
        self.dp50LessCount = 0
        self.dp50MoreCount = 0
        self.dp50Wet = 0
        self.dp50Plus20 = 0
        self.minDpList = []
        self.minDpAsymList = []
        self.maxDpAsymList = []
        self.bList = []
        self.dList = []
        self.cList = []
        self.dp50List = []
        self.dp50WetList = []
        self.dp50Plus20List = []
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
        self.kappaExcel = None
        self.shiftList = []
        self.kappaCalculatedDict = {}
        self.alphaPineneDict = {}
        self.peakPositionInData = []
        self.usableForKappaCalList = []
        self.kappaExcludeList = []
        self.kappaPoints = []
        gc.collect()

    def makeProgress(self, message=None, maxValue=None, complete=False, value=1):
        self.view.makeProgress(message, maxValue, complete, value)
        if self.cancel == True:
            self.cancel = False
            raise InterruptError()

    def cancelProgress(self):
        self.cancel = True

def main():
    #init the controller
    controller = Controller(False)
    #init the view
    view = MainWindow(controller)
    controller.view = view
    view.show_UI()

if __name__ == '__main__':
    main()

=======
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
        self.current_run = 0
        self.number_of_run = 0
        self.run_start_time_list = None
        self.run_end_time_list = None
        self.experiment_date = None
        self.time_of_each_run = 0
        self.ccnc_data = None
        self.smps_data = None
        self.diameter_list = None
        self.dnlog_list = None
        self.ss_list = []
        self.processed_data = None
        self.raw_data = None
        self.current_phase = 0
        self.current_point = None
        self.finish_optimized_phase = False
        self.cancelling = False

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

        self.adjusted_graph = None
        self.dry_diameter_graph = None
        self.min_compare_grpah = None
        self.kappa_graph = None
        self.temp_graph = None
        self.adjusted_graph_list = []
        self.dry_diameter_graph_list = []
        self.min_pos_SMPS_list = []
        self.min_pos_CCNC_list = []
        self.temp_graph_list = []

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
        self.shift_list = []
        self.kappa_exclude_list = []
        self.kappa_points = []
        self.klines = None
        self.max_kappa = None
        self.min_kappa = None
        self.is_full_kappa_graph = False
        self.kappa_calculate_dict = {}
        self.alpha_pinene_dict = {}

        self.run_position_in_data = []
        self.usable_for_kappa_cal_list = []

        self.clone = None


    # -------------General processed_data processing------------------

    def get_flow_rate(self):
        """
        Get the flow rate of the processed_data from user input
        """
        self.flow_rate = self.view.get_flow_rate()
        self.flow_rate = 1/(self.flow_rate*1000/60)

    def get_raw_data(self):
        """
        Get the smps and ccnc processed_data from the input files.
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

        self.ccnc_data = process_csv_files(self.ccnc_csv_files)
        self.smps_data = process_text_files(self.smps_txt_files)
        self.experiment_date = self.ccnc_data[0]
        self.run_start_time_list = self.smps_data[0]
        # Get the running time (in seconds) of each run
        for i in range(len(self.smps_data)):
            if ''.join(self.smps_data[i][0].split()).lower() == "scanuptime(s)":
                scan_up_time = self.smps_data[i][1]
                scan_down_time = self.smps_data[i+1][1]
                break
        self.time_of_each_run = int(scan_up_time) + int(scan_down_time)
        self.run_end_time_list = []
        for i in range(len(self.run_start_time_list)):
            self.run_end_time_list.append(datetime.strptime(self.run_start_time_list[i], "%H:%M:%S") + timedelta(
                seconds=self.time_of_each_run))
        #remove the time frames that doesn't have processed_data in the CCNC
        while (True):
            smps_end_time = self.run_end_time_list[-1]
            csv_end_time = datetime.strptime(self.ccnc_data[1][len(self.ccnc_data[1]) - 1][0], "%H:%M:%S")
            if smps_end_time > csv_end_time:
                self.run_start_time_list = self.run_start_time_list[:-1]
                self.run_end_time_list = self.run_end_time_list[:-1]
            else:
                break
        self.run_end_time_list = [datetime.strftime(x, "%H:%M:%S") for x in self.run_end_time_list]
        self.number_of_run = len(self.run_start_time_list)

    def get_dnlog_data(self):
        """
        Get the dnlog processed_data from the txt file
        """
        for k in range(1, len(self.smps_data)):
            if re.search('[aParam-zA-Z]', self.smps_data[k][0]):
                dnlogPos = k
                break
        startTime = ['dp'] + self.run_start_time_list
        endTime = ['dp'] + self.run_end_time_list
        dnlog_list = [startTime, endTime]
        for k in range(1, dnlogPos):
            dnlog_list.append(self.smps_data[k][:len(startTime)])
        self.dnlog_list = dnlog_list

    def get_smps_and_ccnc_data(self):
        """
        Get the SMPS processed_data from raw SMPS and CCN from raw CCN
        """
        # Get the processed_data from the SMPS file
        ccnc_data = self.ccnc_data[1]
        smps_data = self.smps_data
        startTime = self.run_start_time_list
        endTime = self.run_end_time_list
        width = self.number_of_run

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
        self.run_position_in_data = []
        while (True):
            sizeSum = 0
            countSum = 0
            previousTimeStamp = datetime.strptime(startTime[i], "%H:%M:%S").time()
            maxSS = 0
            # Get all processed_data within the time frame
            self.run_position_in_data.append(len(smpsCcnList))
            for t in range(self.time_of_each_run):
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
            k += self.time_of_each_run

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
                self.run_start_time_list = startTime
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
        """
        Calculate the stretch factor of the ccnc measurement
        """
        self.usable_for_kappa_cal_list = []
        self.min_pos_CCNC_list = []
        self.min_pos_SMPS_list = []
        additionalDataCount = int(self.time_of_each_run / 2)
        startTime = 0
        endTime = self.time_of_each_run
        currPeak = 0
        newData = [self.raw_data.columns.values.tolist()]
        minDist = self.time_of_each_run / 10
        shiftFactor = 0
        self.shift_list = []
        self.number_of_run = 0
        while True:
            # count current_run smps
            aPeak = numpy.asarray(self.raw_data[startTime:endTime]["SMPS Count"])
            minPosSMPS = get_min_index(aPeak)
            # assuming that count current_run of smps is always correct, the first current_run is in the right position
            if minPosSMPS == -1:
                self.min_pos_SMPS_list.append(None)
                self.min_pos_CCNC_list.append(None)
                self.usable_for_kappa_cal_list.append(False)
                minPosSMPS = 0
                minPosCCNC = 0
            else:
                self.min_pos_SMPS_list.append(minPosSMPS + startTime)
                aPeak = numpy.asarray(self.raw_data[startTime + shiftFactor:endTime + additionalDataCount]["CCNC Count"])
                minPosCCNC = get_min_index(aPeak)
                if minPosCCNC == -1:
                    self.min_pos_SMPS_list[-1] = None
                    self.min_pos_CCNC_list.append(None)
                    self.usable_for_kappa_cal_list.append(False)
                    minPosSMPS = 0
                    minPosCCNC = 0
                else:
                    self.usable_for_kappa_cal_list.append(True)
                    shiftFactor = shiftFactor + minPosCCNC - minPosSMPS
                    if currPeak == 0:
                        self.min_pos_CCNC_list.append(minPosSMPS)
                    else:
                        self.min_pos_CCNC_list.append(minPosCCNC + minPosCCNC - minPosSMPS + startTime)
            self.shift_list.append(shiftFactor)

            # Add all peaks to processed_data, whether a peak is a good one or not.
            for i in range(self.time_of_each_run):
                scanTime = i + 1
                realTime = self.raw_data.iat[startTime + i, 1]
                dp = self.raw_data.iat[startTime + i, 2]
                smpsCount = self.raw_data.iat[startTime + i, 3]
                ccncCount = self.raw_data.iat[startTime + i + shiftFactor, 4]
                aveSize = self.raw_data.iat[startTime + i + shiftFactor, 5]
                t1 = self.raw_data.iat[startTime + i + shiftFactor, 6]
                t2 = self.raw_data.iat[startTime + i + shiftFactor, 7]
                t3 = self.raw_data.iat[startTime + i + shiftFactor, 8]
                newData.append([scanTime, realTime, dp, smpsCount, ccncCount, aveSize,t1,t2,t3])

            currPeak += 1
            self.number_of_run += 1
            if endTime + self.time_of_each_run + shiftFactor >= len(self.raw_data) or currPeak >= len(self.run_position_in_data):
                break

            startTime = self.run_position_in_data[currPeak]
            endTime = startTime + self.time_of_each_run
        headers = newData.pop(0)
        self.processed_data = pandas.DataFrame(newData, columns=headers)
        self.processed_data = convert_date(self.processed_data)

    def create_run_alignment_graph(self):
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
                    tempSMPSPeakList.append(tempSMPSPeakList[-1] + self.time_of_each_run)
                    tempCCNPeakCList.append(tempCCNPeakCList[-1] + self.time_of_each_run)
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
        if peak != self.current_run:
            self.switchToPeak(peak)

    # -------------Single peak processing ------------------

    def prepare_run_data(self):
        """
        Convert processed_data to separate lists for easier access. Also normalize and calculate necessary lists.
        """
        try:
            ccnCNList = []
            dNdLogDpList = []
            startPoint = self.current_run * self.time_of_each_run
            endPoint = (self.current_run + 1) * self.time_of_each_run
            self.diameter_list = list(self.processed_data[startPoint:endPoint]['dp'])
            ccnList = list(self.processed_data[startPoint:endPoint]['CCNC Count'])
            cnList = list(self.processed_data[startPoint:endPoint]['SMPS Count'])
            self.temp1 = list(self.processed_data[startPoint:endPoint]['T1'])
            self.temp2 = list(self.processed_data[startPoint:endPoint]['T2'])
            self.temp3 = list(self.processed_data[startPoint:endPoint]['T3'])

            # modify this to get aParam decently good number
            cnList = [x * 0.2 for x in cnList]
            self.ccn_list = ccnList
            self.cn_list = cnList
            self.diameter_midpoint_list = []
            self.ccncnList = ccnCNList
            for i in range(len(ccnList)):
                self.ccncnList.append(ccnList[i] / cnList[i])
            self.dropSizeList = list(self.processed_data[startPoint:endPoint]['Ave Size'])
            for i in range(2, len(self.dnlog_list)):
                self.diameter_midpoint_list.append(self.dnlog_list[i][0])
                dNdLogDpList.append(self.dnlog_list[i][self.current_run + 1])
            self.ccn_normalized_list = normalize_list(dNdLogDpList)
            self.usable_for_kappa_cal_list[self.current_run] = True
        except:
            raise OptimizationError()

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
        """
        Perform charge correction
        """
        try:
            # If the peak doesn't have good alignment, do nothing
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
            frac1List = calculate_fraction(self.diameter_list, 1, coeficientList[1])
            frac2List = calculate_fraction(self.diameter_list, 2, coeficientList[2])
            frac3List = calculate_fraction(self.diameter_list, 3)
            chargeList = []

            for i in self.diameter_list:
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
            maxUpperBinBound = (self.diameter_list[-1] + self.diameter_list[-2]) / 2
            lenDpList = len(self.diameter_list)

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
                        upperBinBound = (self.diameter_list[j] + self.diameter_list[j + 1]) / 2
                        lowerBinBound = (self.diameter_list[j] + self.diameter_list[j - 1]) / 2
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
                        upperBinBound = (self.diameter_list[j] + self.diameter_list[j + 1]) / 2
                        lowerBinBound = (self.diameter_list[j] + self.diameter_list[j - 1]) / 2
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
            self.usable_for_kappa_cal_list[self.current_run] = True
        except:
            raise OptimizationError()

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
                asymList = get_asym_list(self.diameter_list, self.ccnc_sig_list)
                mpIndex = 0

                # Find the index of midPoint
                checkLength = self.time_of_each_run / 20
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
                        self.min_dp = self.diameter_list[i]
                        break

                # Get minDpAsym
                for i in range(mpIndex, mpIndex + self.time_of_each_run / 10):
                    if i < len(self.diameter_list):
                        if self.ccnc_sig_list[i] > currMax:
                            minDpAsymPos = i
                            currMax = self.ccnc_sig_list[i]
                            self.min_dp_asym = self.diameter_list[i]

                maxCCN = max(self.ccnc_sig_list)
                if maxCCN <= 1.2:
                    stableThreshold = 0.075
                else:
                    stableThreshold = 0.25
                # determine maxDplen(self.ccncSigList)
                for i in range(minDpAsymPos + 5, len(self.ccnc_sig_list)):
                    if self.ccnc_sig_list[i] > 1.3:
                        self.max_dp_asym = self.diameter_list[i]
                        break
                    elif abs(self.ccnc_sig_list[i] - self.ccnc_sig_list[i - 1]) > 2 * stableThreshold:
                        self.max_dp_asym = self.diameter_list[i - 1]
                        break

            self.max_dp = self.max_dp_asym
            # Get the processed_data
            asymsList = []
            for i in range(len(self.diameter_list)):
                if self.min_dp_asym < self.diameter_list[i] < self.max_dp_asym:
                    asymsList.append(self.ccnc_sig_list[i])
                else:
                    asymsList.append(0)
            # Calcualte constants
            self.b = get_ave_none_zero(asymsList)
            self.ccn_cn_sim_list.append(0)
            for i in range(1, len(self.diameter_list)):
                if self.min_dp < self.diameter_list[i] < self.max_dp:
                    n = self.b / (1 + (self.diameter_list[i] / self.d) ** self.c)
                    self.ccn_cn_sim_list.append(n)
                else:
                    self.ccn_cn_sim_list.append(self.ccn_cn_sim_list[i - 1])
            self.usable_for_kappa_cal_list[self.current_run] = True
        except:
            raise OptimizationError()

    def fit_sigmoid_line_to_run(self):
        """
        fit the sigmoid line to the the data points.
        """
        try:
            xList = []
            yList = []
            for i in range(len(self.diameter_list)):
                if self.min_dp < self.diameter_list[i] < self.max_dp:
                    xList.append(self.diameter_list[i])
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
            for i in range(1, len(self.diameter_list)):
                if self.diameter_list[i] > self.d:
                    self.dp50_wet = self.dropSizeList[i - 1]
                    break
            for i in range(1, len(self.diameter_list)):
                if self.diameter_list[i] > (self.d + 20):
                    self.dp50_plus_20 = self.dropSizeList[i - 1]
                    break

            for i in range(1, len(self.diameter_list)):
                if self.min_dp < self.diameter_list[i] < self.max_dp:
                    n = self.b / (1 + (self.diameter_list[i] / self.d) ** self.c)
                    self.ccn_cn_sim_list.append(n)
                else:
                    self.ccn_cn_sim_list.append(self.ccn_cn_sim_list[i - 1])
            self.usable_for_kappa_cal_list[self.current_run] = True
        except:
            raise OptimizationError()


    def create_adjusted_graph(self, new_figure = None):
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
            x = range(self.time_of_each_run)
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

    def create_ccn_graph(self, new_figure = None):
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
            plt.plot(self.diameter_list, self.ccncnList, 'o', color="#2196F3", mew=0.5,mec = "#0D47A1",  ms = 9, label = "CCN/CN")
            yLim = min(2, max(self.ccncnList)) + 0.2
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

            x = range(self.time_of_each_run)
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
            self.temp_graph_list.append(plt.gcf())

    def new_figure(self, new_figure = None):
        """
        Make complete graph of the dry diameter after optimization and sigmodal fit
        """
        try:
            if not self.min_pos_CCNC_list[self.current_run] or not self.min_pos_SMPS_list[self.current_run]:
                self.create_ccn_graph(new_figure)
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
            # Graph the sig fit only if the peak is usable
            if self.usable_for_kappa_cal_list[self.current_run]:
                plt.plot(self.diameter_list, self.ccn_cn_sim_list, linewidth=5, color='#EF5350', label="Sigmodal Fit")
            plt.plot(self.diameter_list, self.ccncnList, 'o', color="#2196F3", mew=0.5,mec = "#1976D2",  ms = 9, label = "CCN/CN")
            plt.plot(self.diameter_list, self.ccnc_sig_list, 'o', color="#1565C0", mew=0.5, mec ="#0D47A1", ms = 9, label ="CCN/CN (Corrected)")
            plt.xlabel("Dry diameter(nm)")
            plt.ylabel("CCN/CN ratio and Normalized dN/dLogDP")
            handles, labels = plt.gca().get_legend_handles_labels()
            legend = plt.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.7, 1.1))
            legend.get_frame().set_facecolor('#9E9E9E')
            self.dry_diameter_graph_list.append(plt.gcf())
        except:
            figure = plt.figure(facecolor=settings.graphBackgroundColor)
            self.dry_diameter_graph_list.append(plt.gcf())

    def create_init_graphs(self):
        """
        Make initial graph without optimizations
        """
        plt.ioff()
        self.create_adjusted_graph()
        self.create_ccn_graph()
        self.create_temperature_graph()

    def create_optimized_graphs(self):
        """
        Make complete graphs after optimization
        """
        plt.ioff()
        self.new_figure()

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

    # -------------Single peak procedure ------------------

    def process_one_peak(self):
        """
        The steps required to process a single peak
        """
        try:
            self.prepare_run_data()
            self.make_progress()
            self.create_init_graphs()
            self.make_progress()
        except OptimizationError:
            # Make the peak invalid
            self.min_pos_CCNC_list[self.current_run] = None
            self.min_pos_SMPS_list[self.current_run] = None
            self.usable_for_kappa_cal_list[self.current_run] = False

    def optimize_one_run(self):
        """
        optimize a single run
        """
        try:
            # If a peak is invalid, then do nothing
            if not self.min_pos_CCNC_list[self.current_run] or not self.min_pos_SMPS_list[self.current_run]:
                self.prepare_run_data()
                raise OptimizationError()
            else:
                self.ccnc_sig_list = []
                self.prepare_run_data()
                self.make_progress()
                remove_small_ccn(self.ccn_list, self.min_ccn)
                self.make_progress()
                self.init_correct_charges()
                for i in range(5):
                    self.correct_charges()
                    self.make_progress()
                for i in range(len(self.ccn_fixed_list)):
                    self.ccnc_sig_list.append(self.ccn_fixed_list[i] / self.cn_fixed_list[i])
                self.make_progress()
                self.get_constants()
                self.make_progress()
                self.fit_sigmoid_line_to_run()
                self.make_progress()
                self.create_optimized_graphs()
                self.make_progress()
        except OptimizationError:
            # Disable the peak
            self.usable_for_kappa_cal_list[self.current_run] = False
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
            self.new_figure()
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
            self.dp50_list.append((self.d, self.ss_list[self.current_run]))

    # --------------General processing procedure-------------

    def preprocess_data(self):
        """
        Pre-process the processed_data. This includes extract the dnlog processed_data, smps and ccnc processed_data
        from the raw processed_data file.
        """
        try:
            # reset all constants and the GUI
            self.reset()
            self.view.reset()
            #initiate the progress bar
            self.make_progress(max_value=5)
            self.get_flow_rate()
            self.make_progress("Processing raw processed_data from files...")
            self.get_raw_data()
            self.make_progress("Acquiring dnlog processed_data...")
            self.get_dnlog_data()
            self.make_progress("Acquiring SMPS and CCNC processed_data...")
            self.get_smps_and_ccnc_data()
            self.make_progress("Transforming the CCNC processed_data to match SMPS processed_data....")
            self.match_smps_ccnc_data()
            self.current_phase = 1
            self.make_progress(complete=True)
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

    def initialProcessProcedure(self):
        """
        The processes required to process all the peak processed_data
        """
        try:
            if self.current_phase < 1:
                self.view.show_error_dialog("Data preparation process is not completed.")
                raise InterruptError
            self.make_progress("Preparing processed_data for peak processing...", max_value=self.number_of_run * 3)
            for i in range(0, self.number_of_run):
                self.current_run = i
                self.make_progress("Processing peak " + str(i + 1), value=0)
                self.process_one_peak()
            self.current_run = -1
            self.make_progress(complete=True)
            self.view.update_general_information()
            self.create_run_alignment_graph()
            self.switchToPeak(0)
            self.current_phase = 2
        except InterruptError:
            self.view.show_error_dialog("The processed_data initialization process is cancelled")

    def optimizationProcedure(self):
        """
        The processes required to optimize/sigmodal fit all peaks
        """
        try:
            if self.current_phase < 2:
                self.view.show_error_dialog("Data initialization process is not completed.")
                return

            #clear temp graphs
            for aFigure in self.temp_graph_list:
                aFigure.clf()
                plt.close(aFigure)
            self.temp_graph_list = []
            # clear dry diameter graphs
            for aFigure in self.dry_diameter_graph_list:
                aFigure.clf()
                plt.close(aFigure)
            self.dry_diameter_graph_list = []
            self.make_progress("Preparing processed_data for peak processing...", max_value=self.number_of_run * 12)
            self.finish_optimized_phase = True
            for i in range(0, self.number_of_run):
                self.current_run = i
                self.make_progress("Processing peak " + str(i + 1), value=0)
                self.optimize_one_run()
            # reset current_run
            self.current_run = -1
            self.current_phase = 3
            self.make_progress(complete=True)
            self.switchToPeak(0)
        except InterruptError:
            self.view.show_error_dialog("The optimization process is cancelled")
            self.finish_optimized_phase = False
            for aFigure in self.temp_graph_list:
                aFigure.clf()
                plt.close(aFigure)
            self.temp_graph_list = []
            # clear dry diameter graphs
            for aFigure in self.dry_diameter_graph_list:
                aFigure.clf()
                plt.close(aFigure)
            self.dry_diameter_graph_list = []
            self.initialProcessProcedure()



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
        self.make_progress("Calculating Kappa...", max_value=len(self.dp50_list) + 3)
        self.make_progress("Reading in lookup table")
        if self.kappa_excel is None:
            self.kappa_excel = pandas.read_csv("kCal.csv", header=None)
        lookup = self.kappa_excel
        self.make_progress("Calculating basic consts...")
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
                self.make_progress("Calculating sc...")
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
            self.make_progress("Processing peak" + str(i + 1))
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
        self.make_progress(complete=True)

    #--------------- UI interaction - Data interaction functions--------------------

    def updateGUI(self):
        """
        Update the GUI with the newest processed_data of the current peak
        """
        self.adjusted_graph = self.adjusted_graph_list[self.current_run]
        self.dry_diameter_graph = self.dry_diameter_graph_list[self.current_run]
        self.view.update_dp_dnlog_figures(self.adjusted_graph, self.dry_diameter_graph)
        if self.finish_optimized_phase:
            self.b = self.b_list[self.current_run]
            self.d = self.d_list[self.current_run]
            self.c = self.c_list[self.current_run]
            self.min_dp_asym = self.min_dp_asym_list[self.current_run]
            self.min_dp = self.min_dp_list[self.current_run]
            self.max_dp_asym = self.max_dp_asym_list[self.current_run]
            self.super_saturation_rate = self.dp50_list[self.current_run][1]
            self.dp50 = self.dp50_list[self.current_run][0]
            self.dp50_wet = self.dp50_wet_list[self.current_run]
            self.dp50_plus_20 = self.dp50_plus_20_list[self.current_run]
            self.dp50_less_count = 0
            self.dp50_more_count = 0
            for i in self.diameter_list:
                if i >= self.dp50:
                    self.dp50_more_count += 1
                else:
                    self.dp50_less_count += 1
            self.current_point.set_xdata(numpy.asarray(self.min_pos_SMPS_list)[self.current_run])
            self.current_point.set_ydata(numpy.asarray(self.min_pos_CCNC_list)[self.current_run])
            self.view.update_temp_and_min_figure(self.min_compare_grpah)
            self.view.update_information_after_sig_fit()

        else:
            self.temp_graph = self.temp_graph_list[self.current_run]
            self.super_saturation_rate = self.ss_list[self.current_run]
            self.view.update_temp_and_min_figure(self.temp_graph)
            self.view.update_information_of_run()

    def shiftOneSecond(self, forward=True):
        try:
            if self.current_phase <=1:
                return
            # If the peak is invalid, do nothing
            if not self.usable_for_kappa_cal_list[self.current_run] or self.min_pos_SMPS_list[self.current_run] is None or \
                            self.min_pos_CCNC_list[self.current_run] is None:
                return
            # Get the shift factor
            shiftFactor = self.shift_list[self.current_run]
            # Change the shift factor
            if forward == True:
                shiftFactor -= 1
            else:
                shiftFactor += 1
            self.shift_list[self.current_run] = shiftFactor
            # Get the new processed_data list from raw processed_data
            startTime = self.run_position_in_data[self.current_run] + shiftFactor
            endTime = startTime + self.time_of_each_run
            newCCN = list(self.raw_data.iloc[startTime:endTime, 4])
            newAveSize = list(self.raw_data.iloc[startTime:endTime, 5])
            newT1 = list(self.raw_data.iloc[startTime:endTime, 6])
            newT2 = list(self.raw_data.iloc[startTime:endTime, 7])
            newT3 = list(self.raw_data.iloc[startTime:endTime, 8])

            # Get the location of the replacement in the processed_data
            dataStartTime = self.time_of_each_run * self.current_run
            dataEndTime = dataStartTime + self.time_of_each_run

            # Replace the processed_data in the processed_data
            for i in range(self.time_of_each_run):
                self.processed_data.iat[dataStartTime + i, 4] = newCCN[i]
                self.processed_data.iat[dataStartTime + i, 5] = newAveSize[i]
                self.processed_data.iat[dataStartTime + i, 6] = newT1[i]
                self.processed_data.iat[dataStartTime + i, 7] = newT2[i]
                self.processed_data.iat[dataStartTime + i, 8] = newT3[i]

            # Update the processed_data
            self.prepare_run_data()
            # Fix the graph
            figure = self.adjusted_graph
            self.create_adjusted_graph(new_figure=figure)
            figure = self.dry_diameter_graph
            self.create_ccn_graph(figure)
            figure = self.temp_graph
            self.create_temperature_graph(figure)
            # Update the graph in view
            self.view.updateDpdnlogFigures(self.adjusted_graph, self.dry_diameter_graph)
            self.view.update_temp_and_min_figure(self.temp_graph)
        except:
            self.view.show_error_dialog("Can't shift processed_data. You should disable this peak!")

    def disablePeak(self):
        # disable the peak
        if self.current_phase == 0:
            return
        if not self.finish_optimized_phase:
            self.min_pos_SMPS_list[self.current_run] = None
            self.min_pos_CCNC_list[self.current_run] = None
            self.usable_for_kappa_cal_list[self.current_run] = False
        else:
            self.usable_for_kappa_cal_list[self.current_run] = False
        self.updateGUI()

    def switchToPeak(self, peak):
        """
        Switch to a new peak
        :param peak: the new peak
        """
        if peak != self.current_run:
            self.current_run = peak
            self.updateGUI()

    def reOptimization(self, minDp, minDpAsym, maxDpAsym):
        """
        Re-optimize the processed_data after changing constants minDp, minDpAsym, maxDpAsym
        """
        if not self.finish_optimized_phase:
            return
        self.make_progress("Re-Optimizaing peak" + str(self.current_run + 1), max_value=6)
        self.usable_for_kappa_cal_list[self.current_run] = True
        self.ccnc_sig_list = []
        # Get processed_data
        try:
            self.prepare_run_data()
            self.make_progress()
            remove_small_ccn(self.ccn_list, self.min_ccn)
            self.make_progress()
            self.init_correct_charges()
            for i in range(5):
                self.correct_charges()
                self.make_progress()
            for i in range(len(self.ccn_fixed_list)):
                self.ccnc_sig_list.append(self.ccn_fixed_list[i] / self.cn_fixed_list[i])
            self.make_progress()
            # Update the new constants and reoptimize
            self.get_constants(minDp, minDpAsym, maxDpAsym)
            self.make_progress()
            self.fit_sigmoid_line_to_run()
            self.make_progress()
            # Remake the graph
            figure = self.dry_diameter_graph
            self.new_figure(figure)

            # Update the old processed_data
            self.min_dp_list[self.current_run] = self.min_dp
            self.min_dp_asym_list[self.current_run] = self.min_dp_asym
            self.max_dp_asym_list[self.current_run] = self.max_dp_asym
            self.b_list[self.current_run] = self.b
            self.d_list[self.current_run] = self.d
            self.c_list[self.current_run] = self.c
            self.dp50_list[self.current_run] = (self.d, self.ss_list[self.current_run])
            self.updateGUI()
            self.make_progress(complete=True)
        except:
            # Disable the peak
            self.usable_for_kappa_cal_list[self.current_run] = False
            self.min_dp_list[self.current_run] = 0
            self.min_dp_asym_list[self.current_run] = 0
            self.max_dp_asym_list[self.current_run] = 0
            self.b_list[self.current_run] = 0
            self.d_list[self.current_run] = 0
            self.c_list[self.current_run] = 0
            self.dp50_list[self.current_run] = (0, 0)
            self.updateGUI()
            self.make_progress(complete=True)

    #----------------------Misc------------------------------------------

    def run(self):
        """start to process the processed_data"""
        self.current_phase = 0
        try:
            self.preprocess_data()
        except DataPreparationError:
            pass
        else:
            self.initialProcessProcedure()

    def reset(self):
        """
        Reset all values of the program to process a new processed_data set
        Clear out the memory as well
        """
        i = 1
        gc.collect()
        for aFigure in self.temp_graph_list:
            aFigure.clf()
            plt.close(aFigure)
        self.temp_graph_list = []
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
        self.finish_optimized_phase = False
        self.cancelling = False
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
        self.shift_list = []
        self.kappa_calculate_dict = {}
        self.alpha_pinene_dict = {}
        self.run_position_in_data = []
        self.usable_for_kappa_cal_list = []
        self.kappa_exclude_list = []
        self.kappa_points = []
        gc.collect()

    def make_progress(self, message=None, max_value=None, complete=False, value=1):
        self.view.make_progress(message, max_value, complete, value)
        if self.cancelling == True:
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

>>>>>>> 3cba808e982eb6782e435123a14c81d3423ea1d5
