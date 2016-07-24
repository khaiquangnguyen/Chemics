import csv
import re
import pandas
from ggplot import *
from datetime import *
from scipy import signal
import numpy
import peakutils
from scipy import *
from scipy import stats
from GUI import *
import tempfile
import scipy.constants
import matplotlib.pyplot as plt
import matplotlib
from sys import exit
import settings
matplotlib.style.use('ggplot')
import scipy.optimize as opt
import time
from PySide import QtGui
from Exceptions import *
from HelperFunctions import *
from HelperFunctions import csvProcessing

csvFilePath = "E:\Updated program\Demonstration_Files\CCN data 100203092813.csv"
txtFilePath = "E:\Updated program\Demonstration_Files\AS_Calibration_SMPS.txt"
b = 0

class Controller():
    def __init__(self, mode=True):
        """"
        The controller class for the program. ALso handle the task of the model
        """
        # True if testing, false if not
        self.mode = mode
        # The view associated with the controller
        self.view = None
        # The files which contains the data files
        self.files = None
        # The list of paths to the txt files
        self.smpsTxtFilePath = []
        # The list of paths to the cccn files
        self.ccncCsvFilePath = []
        # The current analysing currPeak
        self.currPeak = 0
        # The total number of peaks
        self.maxPeak = 0
        # The single smps file which contains the data
        self.smpsFile = None
        # The single ccnc file which contains the data, created by combining all the ccn files
        self.ccncFile = None
        # The temporary files used to store the data needed for the program to function
        self.tempFolder = tempfile.mkdtemp('Temp')
        # The list of start time entries
        self.startTimeEntries = None
        # the list of end time entries
        self.endTimeEntries = None
        # The date of the experiment
        self.date = None
        # the time frame of aParam aParam measurement
        self.timeFrame = 0
        # The content of the csv file
        self.csvContent = None
        # The content of the txt file
        self.txtContent = None
        # The list of all diameter of measurement
        self.diameterList = None
        # The DNlog data
        self.dNlog = None
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

        # Necessary information
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

        self.ccnNormalizedFullList = []
        self.ccncnFullList = []
        self.ccncnSimFullList = []
        self.ccncSigFullList = []

        # Variables for peak analysis
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

        # Graph variables
        self.adjustedGraph = None
        self.dryDiaGraph = None
        self.minCompareGraph = None
        self.kappaGraph = None
        self.adjustedGraphList = []
        self.dryDiaGraphList = []
        self.minPosSMPSList = []
        self.minPosCCNCList = []

        # Variables for calculating kappa
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
        self.usablePeakList = []

    #-------------File processing------------------

    def selectFiles(self, files):
        """
        Set the files where data is stored
        :param files: path to files
        """
        self.completedStep = 0
        self.files = files
        self.run()

    def getFileNames(self):
        """
        Get all the SMPS and CCNC file names from the directory
        """
        self.smpsTxtFilePath = []
        self.ccncCsvFilePath = []
        for aFile in self.files:
            if aFile.lower().endswith('.txt'):
                self.smpsTxtFilePath.append(aFile)
            elif aFile.lower().endswith('.csv'):
                self.ccncCsvFilePath.append(aFile)

    def mergeCSVFiles(self):
        """
        If there are more than 1 CSV file selected, merges the csv files together
        """
        self.smpsTxtFilePath = [str(x) for x in self.smpsTxtFilePath]
        self.smpsTxtFilePath = self.smpsTxtFilePath[0]
        self.ccncCsvFilePath = [str(x) for x in self.ccncCsvFilePath]
        if len(self.smpsTxtFilePath) < 1:
            raise FileNotFoundError()
        if len(self.ccncCsvFilePath) < 1:
            raise FileNotFoundError()

    # -------------General data processing------------------

    def getRawDataFromFiles(self):
        """
        Acquire raw data from the SMPS and CCNC files
        """
        # Get the basic contents
        try:
            self.csvContent = csvProcessing(self.ccncCsvFilePath)
            self.date = self.csvContent[0]
            self.txtContent = txtProcessing(self.smpsTxtFilePath)
            self.startTimeEntries = self.txtContent[0]
            self.endTimeEntries = []

            # Get the running time( in seconds) of each run
            for i in range(len(self.txtContent)):
                if ''.join(self.txtContent[i][0].split()).lower() == "scanuptime(s)":
                    scanUpTime = self.txtContent[i][1]
                    scanDownTime = self.txtContent[i+1][1]
            self.timeFrame = int(scanUpTime) + int(scanDownTime)

            # Get the time stamp of each run
            for i in range(len(self.startTimeEntries)):
                self.endTimeEntries.append(datetime.strptime(self.startTimeEntries[i], "%I:%M:%S") + timedelta(
                seconds=self.timeFrame))
            smpsEndTime = self.endTimeEntries[-1]
            csvEndTime = datetime.strptime(self.csvContent[1][len(self.csvContent[1]) - 1][0], "%I:%M:%S")

            # If there are more data in the SMPS file than in the CCNC file. Remove the last smps data
            if smpsEndTime > csvEndTime:
                self.startTimeEntries = self.startTimeEntries[:-1]
                self.endTimeEntries = self.endTimeEntries[:-1]
            self.endTimeEntries = [datetime.strftime(x,"%I:%M:%S") for x in self.endTimeEntries]

            # Get the total number of run
            self.maxPeak = len(self.startTimeEntries)
        except:
            raise FileProcessingError()

    def getDNlog(self):
        """
        Get the DNlog data from the txt file
        """
        try:
            for k in range(1, len(self.txtContent)):
                if re.search('[aParam-zA-Z]', self.txtContent[k][0]):
                    dNlogPos = k
                    break
            startTime = ['dp'] + self.startTimeEntries
            endTime = ['dp'] + self.endTimeEntries
            dNlogList = [startTime, endTime]
            for k in range(1, dNlogPos):
                dNlogList.append(self.txtContent[k][:len(startTime)])
            self.dNlog = dNlogList
        except Exception:
            raise DNlogDataError()

    def getSMPSAndCCNC(self):
        """
        Get the SMPS data from raw SMPS and CCN from raw CCN
        """
        try:
            # Get the data from the SMPS file
            csvContent = self.csvContent[1]
            txtContent = self.txtContent
            startTime = self.startTimeEntries
            endTime = self.endTimeEntries
            width = self.maxPeak

            # get the count position
            for k in range(3, len(self.txtContent)):
                if re.search('[aParam-zA-Z]', self.txtContent[k][0]):
                    for j in range(k, len(self.txtContent)):
                        if not re.search('[aParam-zA-Z]', self.txtContent[j][0]):
                            countPos = j
                            break
                    break

            # Get the count for each scan smps
            loopCount = 0
            countList = []
            smpsList = []
            sumSMPS = 0
            count = 0

            for k in range(countPos, len(txtContent) - 1):
                loopCount += 1
                for j in range(0, width):
                    sumSMPS += float(txtContent[k][j * 2 + 1]) * int(txtContent[k][j * 2 + 2])
                if not countList:
                    for j in range(1, width + 1):
                        num = int(txtContent[k][j * 2])
                        countList.append(num)
                else:
                    for j in range(1, width + 1):
                        num = int(txtContent[k][j * 2])
                        countList[j - 1] += num

                if loopCount == 10:
                    count += 1
                    smpsList.append([sumSMPS / sum(countList)] + [count] + countList)
                    loopCount = 0
                    sumSMPS = 0
                    countList = []

            # Get the CCNC data
            ccnList = []
            aveSizeList = []
            extraCCNList = []
            extraAveSizeList = []
            previousK = 0
            sizeList = [0.625] + [0.875]
            size = 1.25
            binPos = 25
            aSSList = []

            # Calculate the size of each bin
            for i in range(0, width):
                sizeList.append(size)
                size += 0.5

            # Get the first time stamp position in the ccnc file
            timeStamp = startTime[0]
            for k in range(0, len(csvContent)):
                if csvContent[k][0] == timeStamp:
                    previousK = k
                    break

            i = 0
            t = 0
            isPeakTimeStamp = True
            smpsCcnList = []
            timeStamp = datetime.strptime(startTime[0], "%I:%M:%S")
            self.peakPositionInData = []
            while (True):
                sizeSum = 0
                countSum = 0
                previousTimeStamp = datetime.strptime(startTime[i], "%I:%M:%S").time()
                maxSS = 0
                # Get all data within the time frame
                self.peakPositionInData.append(len(smpsCcnList))
                for t in range(self.timeFrame):
                    aLine = [smpsList[t][1]] + [timeStamp.time()] + [float(smpsList[t][0])] + [
                        float(smpsList[t][i + 2])]
                    sizeSum = 0
                    countSum = 0
                    aLine.append(float(csvContent[k + t][-3]))
                    if csvContent[k + t][1] > maxSS:
                        maxSS = csvContent[k + t][1]
                    for m in range(0, 20):
                        sizeSum += sizeList[m] * float(csvContent[k + t][binPos + m])
                        size += 0.5
                        countSum += float(csvContent[k + t][binPos + m])
                        # get the average size for each scan
                    if countSum == 0:
                        aLine.append(0)
                    else:
                        aLine.append(sizeSum / countSum)
                    smpsCcnList.append(aLine)
                    timeStamp += timedelta(seconds=1)
                k += self.timeFrame
                self.ssList.append(maxSS)

                # If reach the end of smps, collect the rest of ccnc for alignment
                timeGap = 0
                if i == len(startTime) - 1:
                    timeGap = (datetime.strptime(csvContent[-1][0], "%I:%M:%S") - timeStamp).seconds
                # Do whatever
                else:
                    nextTimeStamp = datetime.strptime(startTime[i + 1], "%I:%M:%S")
                    if nextTimeStamp < timeStamp:
                        timeGap = -(timeStamp - nextTimeStamp).seconds
                    else:
                        timeGap = (nextTimeStamp - timeStamp).seconds

                # If the time shown in header is smaller than the supposed actual time frame
                # Then correct the header accordingly
                if timeGap < 0:
                    startTime[i+1] = (datetime.strptime(startTime[i + 1], "%I:%M:%S") + timedelta(seconds = abs(timeGap))).strftime( "%I:%M:%S")
                    self.startTimeEntries = startTime
                    timeGamp = 0
                for t in range(timeGap):
                    aLine = [None] + [timeStamp.time()] + [None] + [None]
                    sizeSum = 0
                    countSum = 0
                    aLine.append(float(csvContent[k + t][-3]))
                    if csvContent[k + t][1] > maxSS:
                        maxSS = csvContent[k + t][1]
                    for m in range(0, 20):
                        sizeSum += sizeList[m] * float(csvContent[k + t][binPos + m])
                        size += 0.5
                        countSum += float(csvContent[k + t][binPos + m])
                        # get the average size for each scan
                    if countSum == 0:
                        aLine.append(0)
                    else:
                        aLine.append(sizeSum / countSum)
                    smpsCcnList.append(aLine)
                    timeStamp += timedelta(seconds=1)
                k += timeGap

                # Loop break condition
                i += 1
                if i >= len(startTime):
                    break

            title = ["scan time"] + ["real time"] + ["dp"] + ["SMPS Count"] + ["CCNC Count"] + ["Ave Size"]
            self.rawData = pandas.DataFrame(smpsCcnList, columns=title)
        except:
            raise DataError()

    def matchSMPSCCNCData(self):
        """
        Calculate the stretch factor of the ccnc measurement
        """
        try:
            self.usablePeakList = []
            self.minPosCCNCList = []
            self.minPosSMPSList = []
            additionalDataCount = int(self.timeFrame / 2)
            startTime = 0
            endTime = self.timeFrame
            currPeak = 0
            newData = [self.rawData.columns.values.tolist()]
            minDist = self.timeFrame / 10
            shiftFactor = 0
            self.shiftList = []

            while True:
                # count currPeak smps
                aPeak = numpy.asarray(self.rawData[startTime:endTime]["SMPS Count"])
                minPosSMPS = getMinIndex(aPeak)
                # assuming that count currPeak of smps is always correct, the first currPeak is in the right position
                if minPosSMPS == -1:
                    self.minPosSMPSList.append(None)
                    self.minPosCCNCList.append(None)
                    self.usablePeakList.append(False)
                    minPosSMPS = 0
                    minPosCCNC = 0
                else:
                    self.minPosSMPSList.append(minPosSMPS + startTime)
                    aPeak = numpy.asarray(self.rawData[startTime + shiftFactor:endTime + additionalDataCount]["CCNC Count"])
                    minPosCCNC = getMinIndex(aPeak)
                    if minPosCCNC == -1:
                        self.minPosSMPSList[-1] = None
                        self.minPosCCNCList.append(None)
                        self.usablePeakList.append(False)
                        minPosSMPS = 0
                        minPosCCNC = 0
                    else:
                        self.usablePeakList.append(True)
                        shiftFactor = shiftFactor + minPosCCNC - minPosSMPS
                        if currPeak == 0:
                            self.minPosCCNCList.append(minPosSMPS)
                        else:
                            self.minPosCCNCList.append(minPosCCNC + minPosCCNC - minPosSMPS + startTime)
                self.shiftList.append(shiftFactor)

                # Add all peaks to data, whether a peak is a good one or not.
                for i in range(self.timeFrame):
                    scanTime = i + 1
                    realTime = self.rawData.iat[startTime + i, 1]
                    dp = self.rawData.iat[startTime + i, 2]
                    smpsCount = self.rawData.iat[startTime + i, 3]
                    ccncCount = self.rawData.iat[startTime + i + shiftFactor, 4]
                    aveSize = self.rawData.iat[startTime + i + shiftFactor, 5]
                    newData.append([scanTime, realTime, dp, smpsCount, ccncCount, aveSize])

                currPeak += 1
                if endTime + shiftFactor >= len(self.rawData) or currPeak >= len(self.peakPositionInData):
                    break

                startTime = self.peakPositionInData[currPeak]
                endTime = startTime + self.timeFrame

            headers = newData.pop(0)
            self.data = pandas.DataFrame(newData, columns=headers)
            self.data = dateConvert(self.data)
        except:
            raise DataMatchingError()

    def peakAlignAndGraph(self):
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
        tempSMPSPeak = []
        tempCCNCList = []
        # Plot only valid points
        for i in range(len(self.minPosCCNCList)):
            if self.minPosCCNCList[i] and self.minPosSMPSList[i]:
                tempSMPSPeak.append(self.minPosSMPSList[i])
                tempCCNCList.append(self.minPosCCNCList[i])

        x = numpy.asarray(tempSMPSPeak)
        y = numpy.asarray(tempCCNCList)

        result = scipy.stats.linregress(x, y)
        slope = result[0]
        yIntercept = result[1]

        # Recalculate the position of the smps
        correctedIndexList = []
        for i in range(len(self.data)):
            correctedIndexList.append(round((i * slope + yIntercept) * 10))
        plt.plot(x, x * slope + yIntercept, linewidth = 4, color = '#43A047', label = "Regression line")
        plt.plot(x, y, "o", ms = 10, color = "#43A047",picker=15, mew = 0, label = "Minimum")
        textToShow = str(slope) + "* x" + " + " + str(yIntercept)
        self.currentPoint, = plt.plot(x[0],y[0],'o', color = "#81C784", ms = 12, mew = 0)
        plt.text(x[4], y[3], textToShow, color = "#81C784" )
        plt.xlabel("SMPS minumum point")
        plt.ylabel("CCNC minimum point")

        handles, labels = plt.gca().get_legend_handles_labels()
        legend = plt.legend(handles, labels,loc="upper right", bbox_to_anchor=(1,0.7))
        legend.get_frame().set_facecolor('#9E9E9E')

        if not self.minCompareGraph:
            self.minCompareGraph = plt.gcf()
        else:
            self.minCompareGraph.clf()
            plt.close(self.minCompareGraph)
            self.minCompareGraph = plt.gcf()

    def onPick(self, event):
        """
        When aParam dot on the graph is selected
        """
        peak = event.ind
        if peak != self.currPeak:
            self.updatePeak(peak)

    # -------------Single peak processing ------------------

    def preparePeakData(self):
        """
        Convert data to separate lists for easier access. Also normalize and calculate necessary lists. Prepare data
        for later usage
        """
        try:
            ccnCNList = []
            dNdLogDpList = []
            startPoint = self.currPeak * self.timeFrame
            endPoint = (self.currPeak + 1) * self.timeFrame
            self.diameterList = list(self.data[startPoint:endPoint]['dp'])
            ccnList = list(self.data[startPoint:endPoint]['CCNC Count'])
            cnList = list(self.data[startPoint:endPoint]['SMPS Count'])

            # modify this to get aParam decently good number
            cnList = [x * 0.2 for x in cnList]
            self.ccnList = ccnList
            self.cnList = cnList
            self.diameterMidpointList = []
            self.ccncnList = ccnCNList
            for i in range(len(ccnList)):
                self.ccncnList.append(ccnList[i] / cnList[i])
            self.dropSizeList = list(self.data[startPoint:endPoint]['Ave Size'])
            for i in range(2, len(self.dNlog)):
                self.diameterMidpointList.append(self.dNlog[i][0])
                dNdLogDpList.append(self.dNlog[i][self.currPeak + 1])
            self.ccnNormalizedList = normalizeList(dNdLogDpList)
            self.usablePeakList[self.currPeak] = True
        except:
            # Make the peak invalid
            self.usablePeakList[self.currPeak] = False

    def initCorrectCharges(self):
        """
        Initiate the correct charge procedure
        """
        # If the peak doesn't have good alignment, do nothing
        if self.minPosCCNCList[self.currPeak] and self.minPosSMPSList[self.currPeak]:
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
            if self.minPosSMPSList[self.currPeak] is None or self.minPosCCNCList[self.currPeak] is None:
                return
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
            self.usablePeakList[self.currPeak] = True
        except:
            self.usablePeakList[self.currPeak] = False

    def getConstants(self, minDp = None, minDPAsym = None, maxDpAsym = None):
        """
        Acquire the necessary constants from the data to perform sigmodal fit
        """
        try:
            # If nothing to get constants from
            if self.minPosCCNCList[self.currPeak] is None or self.minPosSMPSList[self.currPeak] is None:
                return

            if minDp and minDPAsym and maxDpAsym:
                self.minDp = minDp
                self.minDpAsym = minDPAsym
                self.maxDpAsym = maxDpAsym
            else:
                asymList = getAsym(self.diameterList, self.ccncSigList)
                mpIndex = 0

                # Find the index of midPoint
                checkLength = self.timeFrame / 20
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
                for i in range(mpIndex, mpIndex + self.timeFrame / 10):
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
                    elif abs(self.ccncSigList[i] - self.ccncSigList[i-1]) > 2 * stableThreshold:
                        self.maxDpAsym = self.diameterList[i-1]
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

            self.usablePeakList[self.currPeak] = True
        except:
            self.usablePeakList[self.currPeak] = False

    def optimizePeak(self):
        """
        Perform optimization, which is basically sigmodal fit
        """
        try:
            if not self.minPosCCNCList[self.currPeak] or not self.minPosSMPSList[self.currPeak]:
                self.b = 0
                self.d = 0
                self.c = 0
                if len(self.bList) >= self.currPeak:
                    self.bList.append(0)
                    self.dList.append(0)
                    self.cList.append(0)
                    self.dp50List.append((0, 0))
                else:
                    self.bList[self.currPeak] = 0
                    self.dList[self.currPeak] = 0
                    self.cList[self.currPeak] = 0
                    self.dp50List[self.currPeak]= (0,0)
                return

            xList = []
            yList = []
            for i in range(len(self.diameterList)):
                if self.minDp < self.diameterList[i] < self.maxDp:
                    xList.append(self.diameterList[i])
                    yList.append(self.ccncSigList[i])
            initGuess = [30, -10]
            xList = numpy.asarray(xList)
            yList = numpy.asarray(yList)
            result = opt.curve_fit(f, xList, yList, bounds=([self.minDp, -200], [self.maxDp, -1]), method="trf")
            self.d = result[0][0]
            self.c = result[0][1]
            self.ccncnSimList = [0]
            for i in range(1, len(self.diameterList)):
                if self.minDp < self.diameterList[i] < self.maxDp:
                    n = self.b / (1 + (self.diameterList[i] / self.d) ** self.c)
                    self.ccncnSimList.append(n)
                else:
                    self.ccncnSimList.append(self.ccncnSimList[i - 1])
            self.usablePeakList[self.currPeak] = True
        except:
            self.usablePeakList[self.currPeak] = False

    # -------------Produce graphs ------------------

    def makeAdjustedGraph(self, newFigure = None):
        """
        Make the comparable CCN and CN graph after adjustment of minimum peak
        """
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
        x = range(self.timeFrame)
        plt.plot(x, self.cnList, linewidth = 4, color = '#EF5350', label = "CN")
        plt.plot(x, self.ccnList, linewidth = 4, color = '#2196F3',label = "CCN")
        handles, labels = plt.gca().get_legend_handles_labels()
        legend = plt.legend(handles, labels, loc="upper left", bbox_to_anchor=(0, 0.9))
        legend.get_frame().set_facecolor('#9E9E9E')

        plt.xlabel("Scan time(s)")
        plt.ylabel("Concentration(cm3)")
        self.adjustedGraphList.append(plt.gcf())

    def makeCCNGraph(self, newFigure = None):
        """
        Make graph for the CCN data. Used to cross check for peak alignment.
        """
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
        plt.gca().tick_params(axis='x', color='1', which='both', labelcolor="0.6")
        plt.gca().tick_params(axis='y', color='1', which='both', labelcolor="0.6")
        plt.gca().yaxis.label.set_color('0.6')
        plt.gca().xaxis.label.set_color('0.6')
        plt.plot(self.diameterList, self.ccncnList, 'o', color="#2196F3", mew=0.5,mec = "#0D47A1",  ms = 9, label = "CCN/CN")
        plt.gca().axes.set_ylim([-0.1, 1.2])
        handles, labels = plt.gca().get_legend_handles_labels()
        legend = plt.legend(handles, labels,loc="upper right", bbox_to_anchor=(1,0.9))
        legend.get_frame().set_facecolor('#9E9E9E')
        plt.xlabel("Diameter (nm)")
        plt.ylabel("CCN/CN")
        self.dryDiaGraphList.append(plt.gcf())

    def makeFullDryDiameterGraph(self):
        """
        Make complete graph of the dry diameter after optimization and sigmodal fit
        """
        if not self.minPosCCNCList[self.currPeak] or not self.minPosSMPSList[self.currPeak]:
            self.makeCCNGraph()
            return
        figure = plt.figure(facecolor=settings.graphBackgroundColor)
        plt.axes(frameon=False)
        plt.grid(color='0.5')
        plt.axhline(0, color='0.6', linewidth=4)
        plt.axvline(0, color='0.6', linewidth=4)
        plt.gca().tick_params(axis='x', color='1', which='both', labelcolor="0.6")
        plt.gca().tick_params(axis='y', color='1', which='both', labelcolor="0.6")
        plt.gca().yaxis.label.set_color('0.6')
        plt.gca().xaxis.label.set_color('0.6')
        plt.gca().axes.set_ylim([-0.1, 1.2])

        plt.plot(self.diameterMidpointList, self.ccnNormalizedList, linewidth=4, color='#43A047', label="dN/dLogDP")
        # Graph the sig fit only if the peak is usable
        if self.usablePeakList[self.currPeak]:
            plt.plot(self.diameterList, self.ccncnSimList, linewidth=5, color='#EF5350', label="Sigmodal Fit")
        plt.plot(self.diameterList, self.ccncnList, 'o', color="#2196F3", mew=0.5,mec = "#1976D2",  ms = 9, label = "CCN/CN")
        plt.plot(self.diameterList, self.ccncSigList, 'o', color="#1565C0", mew=0.5,mec = "#0D47A1",  ms = 9, label = "CCN/CN (Corrected)")
        plt.xlabel("Dry diameter(nm)")
        plt.ylabel("CCN/CN ratio and Normalized dN/dLogDP")
        handles, labels = plt.gca().get_legend_handles_labels()
        legend = plt.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.7, 1.1))
        legend.get_frame().set_facecolor('#9E9E9E')

        self.dryDiaGraphList.append(plt.gcf())

    def makeInitialGraphs(self):
        """
        Make initial graph without optimizations
        """
        plt.ioff()
        self.makeAdjustedGraph()
        self.makeCCNGraph()

    def makeOptimizedGraphs(self):
        """
        Make complete graphs after optimization
        """
        plt.ioff()
        self.makeFullDryDiameterGraph()

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
            print "There is an error on peak"

    def singlePeakOptimizationProcedure(self):
        """
        The steps required to optimizePeak a single peak
        :return:
        """
        try:
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
            if len(self.minDpList) >= self.currPeak:
                self.minDpList.append(self.minDp)
                self.minDpAsymList.append(self.minDpAsym)
                self.maxDpAsymList.append(self.maxDpAsym)
            else:
                self.minDpList[self.currPeak] = self.minDp
                self.minDpAsymList[self.currPeak] = self.minDpAsym
                self.maxDpAsymList[self.currPeak] = self.maxDpAsym

            self.makeProgress()
            self.optimizePeak()
            if len(self.bList) >= self.currPeak:
                self.bList.append(self.b)
                self.dList.append(self.d)
                self.cList.append(self.c)
                self.dp50List.append((self.d, self.ssList[self.currPeak]))
            else:
                self.bList=self.b
                self.dList=self.d
                self.cList=self.c
                self.dp50List=(self.d, self.ssList[self.currPeak])

            if len(self.ccnNormalizedFullList) <= (self.currPeak+1) * self.timeFrame:
                self.ccnNormalizedFullList.append(self.ccnNormalizedList)
                self.ccncnFullList.append(self.ccncnList)
                self.ccncSigFullList.append(self.ccncSigList)
                self.ccncnSimFullList.append(self.ccncnSimList)
            else:
                startTime = self.currPeak * self.timeFrame
                endTime = self.currPeak * self.timeFrame
                self.ccnNormalizedFullList[startTime:endTime] = self.ccnNormalizedList
                self.ccncnFullList[startTime:endTime] = self.ccncnList
                self.ccncSigFullList[startTime:endTime] = self.ccncSigList
                self.ccncnSimFullList[startTime:endTime] = self.ccncnSimList

            self.makeProgress()
            self.makeOptimizedGraphs()
            self.makeProgress()
        except OptimizationError:
            # self.view.showError("The data is not optimizable. No optimal soluation found!")
            print "There is an error on peak"

    # --------------General processing procedure-------------

    def preparationProcedure(self):
        """
        The process necessary to read in and pre-process the data
        """
        try:
            # reset the GUI
            self.view.reset()
            # clear all figures
            for aFigure in self.adjustedGraphList:
                aFigure.clf()
                plt.close(aFigure)
            self.adjustedGraphList = []
            for aFigure in self.dryDiaGraphList:
                aFigure.clf()
                plt.close(aFigure)
            self.dryDiaGraphList = []
            self.makeProgress(maxValue=7)
            self.makeProgress("Searching files for CCNC and SMPS files...")
            self.getFileNames()
            self.makeProgress("Merging SMPS and CCNC files...")
            self.mergeCSVFiles()
            self.makeProgress("Processing raw data from files...")
            self.getRawDataFromFiles()
            self.makeProgress("Acquiring DnLog data...")
            self.getDNlog()
            self.makeProgress("Acquiring SMPS and CCNC data...")
            self.getSMPSAndCCNC()
            self.makeProgress("Transforming the CCNC data to match SMPS data....")
            self.matchSMPSCCNCData()
            self.completedStep = 1
            self.makeProgress(complete=True)
        except FileNotFoundError:
            self.view.showError("Can't find SMPS or CCNC data in files: " + str(self.files))
            raise DataPreparationError()
        except FileProcessingError:
            self.view.showError("Can't process the SMPS or CCNC files!")
            raise DataPreparationError()
        except DNlogDataError:
            self.view.showError("Can't process DNlog data from the SMPS file!")
            raise DataPreparationError()
        except DataError:
            self.view.showError("Can't process data from SMPS or CCNC raw data!")
            raise DataPreparationError()
        except DataMatchingError:
            self.view.showError("Can't match SMPS and CCNC data!")
            raise DataPreparationError()
        except InterruptError:
            self.view.showError("The preparation process is cancelled")
            raise DataPreparationError()

    def initialProcessProcedure(self):
        # -1 because we kind of making up the last peak for alignment
        # 5 for 5 step for each peak processing
        try:
            if self.completedStep < 1:
                self.view.showError("Data preparation process is not completed.")
                raise InterruptError
            self.makeProgress("Preparing data for peak processing...", maxValue=self.maxPeak * 3)
            for i in range(0, self.maxPeak):
                self.currPeak = i
                self.makeProgress("Processing peak " + str(i + 1), value=0)
                self.singlePeakProcessingProcedure()
            # reset currPeak
            self.currPeak = -1
            self.completedStep = 2
            self.makeProgress(complete=True)
            self.view.updateData()
            self.peakAlignAndGraph()
            self.updatePeak(0)
            self.view.updateTotalViewFigure(self.minCompareGraph)
        except InterruptError:
            self.view.showError("The data initialization process is cancelled")

    def optimizationProcedure(self):
        """
        The optimization process for all peaks
        """
        try:
            if self.completedStep < 2:
                self.view.showError("Data initialization process is not completed.")
                return
            # clear dry diameter graphs
            for aFigure in self.dryDiaGraphList:
                aFigure.clf()
                plt.close(aFigure)
            self.dryDiaGraphList = []
            self.makeProgress("Preparing data for peak processing...", maxValue=self.maxPeak * 12)
            self.optimized = True

            for i in range(0, self.maxPeak):

                self.currPeak = i
                self.makeProgress("Processing peak " + str(i + 1), value=0)
                self.singlePeakOptimizationProcedure()
            # reset currPeak
            self.currPeak = -1
            self.completedStep = 3
            self.makeProgress(complete=True)
            self.updatePeak(0)
        except InterruptError:
            self.view.showError("The optimization process is cancelled")

    #--------------- UI interaction functions--------------------

    def updatePeakData(self):
        """
        Get the information of the current peak
        """
        if self.optimized:
            self.b = self.bList[self.currPeak]
            self.d = self.dList[self.currPeak]
            self.c = self.cList[self.currPeak]
        self.adjustedGraph = self.adjustedGraphList[self.currPeak]
        self.dryDiaGraph = self.dryDiaGraphList[self.currPeak]
        self.superSaturation = self.ssList[self.currPeak]
        self.dp50 = self.d
        self.dp50LessCount = 0
        self.dp50MoreCount = 0
        for i in self.diameterList:
            if i >= self.dp50:
                self.dp50MoreCount += 1
            else:
                self.dp50LessCount += 1
        for i in range(1,len(self.diameterList)):
            if i > self.dp50:
                self.dp50Wet = self.dropSizeList[i-1]
                break
        for i in range(1, len(self.diameterList)):
            if i > (self.dp50 + 20):
                self.dp50Plus20 = self.dropSizeList[i - 1]
                break

        self.view.updateFigures(self.adjustedGraph, self.dryDiaGraph)
        self.view.updateTotalViewFigure(self.minCompareGraph)
        self.currentPoint.set_xdata(numpy.asarray(self.minPosSMPSList)[self.currPeak])
        self.currentPoint.set_ydata(numpy.asarray(self.minPosCCNCList)[self.currPeak])
        self.view.updatePeakInfo()

    def shiftDataCCNC(self, forward=True):
        # If the peak is invalid, do nothing
        if not self.usablePeakList[self.currPeak] or self.minPosSMPSList[self.currPeak] is None or \
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
        endTime = startTime + self.timeFrame
        newCCN = list(self.rawData.iloc[startTime:endTime, 4])
        newAveSize = list(self.rawData.iloc[startTime:endTime, 5])

        # Get the location of the replacement in the data
        dataStartTime = self.timeFrame * self.currPeak
        dataEndTime = dataStartTime + self.timeFrame

        # Replace the data in the data
        for i in range(self.timeFrame):
            self.data.iat[dataStartTime + i, 4] = newCCN[i]
            self.data.iat[dataStartTime + i, 5] = newAveSize[i]

        # Update the data
        self.preparePeakData()
        # Fix the graph
        figure = self.adjustedGraph
        self.makeAdjustedGraph(newFigure=figure)
        figure = self.dryDiaGraph
        self.makeCCNGraph(figure)
        # Update the graph in view
        self.view.updateFigures(self.adjustedGraph, self.dryDiaGraph)

    def makeKappaGraph(self):
        klines = pandas.read_excel("klines.xlsx", header=1)
        header = klines.columns
        diaList = klines[header[1]]
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

        for i in range(2, len(header)):
            y = klines[header[i]]
            plt.loglog(diaList, y, label=str(header[i]), linewidth = 4)

        handles, labels = plt.gca().get_legend_handles_labels()
        legend = plt.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.15, 1))
        legend.get_frame().set_facecolor('#9E9E9E')
        self.kappaGraph = plt.gcf()

    def calKappa(self):
        self.dp50List = [(66.873131326442845, 0.2), (64.706293297900331, 0.2), (66.426791348408827, 0.2), (65.807043010964122, 0.4), (39.029118190703379, 0.4), (41.656041922784382, 0.4), (42.222353379447377, 0.4), (38.860120694533627, 0.4), (38.779984169692248, 0.4), (29.464779084111022, 0.6), (31.946994836267585, 0.6), (32.297643866436054, 0.6), (32.50404169014837, 0.6), (32.495398001104491, 0.6), (122.45185476098608, 0.8), (25.707116797205551, 0.8), (26.295107828742754, 0.8), (26.584143571968784, 0.8)]
        self.kappaExcel = pandas.read_excel("kCal.xlsx", header=None, sheetname=["lookup","apKappa"])
        lookup = self.kappaExcel["lookup"]
        apKappa = self.kappaExcel['apKappa']
        # Calculate constants (mostly)
        # scCalcsC1 = []
        # sc = list(scCalcs.iloc[1:][4])
        # iKappa = 0.00567
        # self.sc = (max(sc) - 1) * 100
        # dd = 280
        # self.asc = (exp(sqrt(4 * self.aParam ** 3 / (27 * iKappa * (dd * 0.000000001) ** 3))) - 1) * 100
        # sc2 = list(scCalcs.iloc[2:][8])
        # self.sc2 = (max(sc2) - 1) * 100
        # trueSC = list(scCalcs.iloc[2:][1])
        # self.trueSC = (max(sc2) - 1) * 100

        # Calculate each kappa
        for i in range(len(self.dp50List)):
            ss = self.dp50List[i][1]
            dp50 = self.dp50List[i][0]
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
            self.anaKappa = (4 * self.aParam ** 3) / (27 * (dp50 * 0.000000001) ** 3 * log(ss / 100 + 1) ** 2)
            kDevi = (self.appKappa-self.anaKappa)/self.appKappa * 100
            if ss in self.kappaCalculatedDict.keys():
                self.kappaCalculatedDict[ss].append([dp50, self.appKappa, self.anaKappa,kDevi])
            else:
                self.kappaCalculatedDict[ss] = ([[dp50, self.appKappa, self.anaKappa, kDevi]])

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
            self.alphaPineneDict[aKey] = (meanDp, stdDp, meanApp, stdApp, meanAna, stdAna, meanDev, devMean)

    def run(self):
        """
        The main running procedure of the program
        """
        # If not testing
        if not self.mode:
            # If view is not exist, initiate the view
            if not self.view:
                view = MainWindow(self)
                self.setView(view)
                view.run()

            # Run the procedure
            else:
                try:
                    self.preparationProcedure()
                except DataPreparationError:
                    pass
                else:
                    self.initialProcessProcedure()


        else:
            self.getRawDataFromFiles(csvFilePath, txtFilePath)
            self.getDNlog()
            self.getSMPSAndCCNC()
            self.matchSMPSCCNCData()
            self.preparePeakData()

    def updatePeak(self, peak):
        if peak != self.currPeak:
            self.currPeak = peak
            self.updatePeakData()

    def setView(self, view):
        """
        Set the view for the controller
        :param view: the view which associate with the controller
        """
        self.view = view

    def showData(self):
        if not self.mode:
            self.view.showData()

    def makeProgress(self, message=None, maxValue=None, complete=False, value=1):
        if not self.mode:
            self.view.makeProgress(message, maxValue, complete, value)
        else:
            print message
        if self.cancel == True:
            self.cancel = False
            raise InterruptError()

    def cancelProgress(self):
        self.cancel = True

    def disablePeak(self):
        # disable the peak
        self.usablePeakList[self.currPeak] = False

    def enablePeak(self):
        """
        Enable the peak
        """
        self.usablePeakList[self.currPeak] = True

    def reOptimization(self,minDp,minDpAsym,maxDpAsym):
        # Enable the peak
        self.makeProgress("Re-Optimizaing peak" + str(self.currPeak + 1), maxValue=6)
        self.enablePeak()
        self.ccncSigList = []
        # Get other data
        self.preparePeakData()
        # Get Sig list data
        startTime = self.currPeak * self.timeFrame
        endTime = startTime + self.timeFrame
        removeSmallCcn(self.ccnList, self.minCcn)
        self.makeProgress()
        self.ccncSigList = self.ccncSigFullList[startTime:endTime]
        self.makeProgress()
        self.getConstants(minDp,minDpAsym,maxDpAsym)

        # Replace the old data
        if len(self.minDpList) >= self.currPeak:
            self.minDpList.append(self.minDp)
            self.minDpAsymList.append(self.minDpAsym)
            self.maxDpAsymList.append(self.maxDpAsym)
        else:
            self.minDpList[self.currPeak] = self.minDp
            self.minDpAsymList[self.currPeak] = self.minDpAsym
            self.maxDpAsymList[self.currPeak] = self.maxDpAsym

        self.makeProgress()
        self.optimizePeak()
        if len(self.bList) >= self.currPeak:
            self.bList.append(self.b)
            self.dList.append(self.d)
            self.cList.append(self.c)
            self.dp50List.append((self.d, self.ssList[self.currPeak]))
        else:
            self.bList = self.b
            self.dList = self.d
            self.cList = self.c
            self.dp50List = (self.d, self.ssList[self.currPeak])

        if len(self.ccnNormalizedFullList) <= (self.currPeak + 1) * self.timeFrame:
            self.ccnNormalizedFullList.append(self.ccnNormalizedList)
            self.ccncnFullList.append(self.ccncnList)
            self.ccncSigFullList.append(self.ccncSigList)
            self.ccncnSimFullList.append(self.ccncnSimList)
        else:
            startTime = self.currPeak * self.timeFrame
            endTime = self.currPeak * self.timeFrame
            self.ccnNormalizedFullList[startTime:endTime] = self.ccnNormalizedList
            self.ccncnFullList[startTime:endTime] = self.ccncnList
            self.ccncSigFullList[startTime:endTime] = self.ccncSigList
            self.ccncnSimFullList[startTime:endTime] = self.ccncnSimList
        self.makeProgress()

        # Remake the graph
        figure = self.dryDiaGraph
        self.makeCCNGraph(figure)
        # Update the graph in view
        self.updatePeakData()


        #

def main():
    controller = Controller(False)
    controller.run()
    # controller.calKappa()
    # test()

if __name__ == '__main__':
    main()
