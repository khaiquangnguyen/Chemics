import csv
import re
import pandas
from ggplot import *
from datetime import *
from scipy import signal
import numpy
import peakutils
from scipy import *
from settings import *
from GUI import *
import tempfile
import scipy.constants
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('ggplot')
import scipy.optimize as opt
import time
from PySide import QtGui

csvFilePath = "E:\Updated program\Demonstration_Files\CCN data 100203092813.csv"
txtFilePath = "E:\Updated program\Demonstration_Files\AS_Calibration_SMPS.txt"
b = 0
class FileNotFoundError(Exception):
    def __init__(self):
        pass

class FileProcessingError(Exception):
    def __init__(self):
        pass

class DNlogDataError(Exception):
    def __init__(self):
        pass

class DataError(Exception):
    def __init__(self):
        pass

class DataMatchingException(Exception):
    def __init__(self):
        pass

class DataPreparationError(Exception):
    def __init__(self):
        pass

class OptimizationError(Exception):
    def __init__(self):
        pass


class Controller():
    def __init__(self,mode = True):
        """"
        The controller class for the program. ALso handle the task of the model
        """
        # True if testing, false if not
        self.mode = mode
        # The view associated with the controller
        self.view = None
        # The folder which contains the data files
        self.folder = None
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
        # The temporary folder used to store the data needed for the program to function
        self.tempFolder = tempfile.mkdtemp('Temp')
        # The list of start time entries
        self.startTimeEntries = None
        # the list of end time entries
        self.endTimeEntries = None
        # The date of the experiment
        self.date = None
        # the time frame of a a measurement
        self.timeFrame = 0
        # The content of the csv file
        self.csvContent = None
        # The content of the txt file
        self.txtContent = None
        # The list of all diameter of measurement
        self.diameterList = None
        # The DNlog data
        self.dNlog = None
        # The SMPS and CCNC data
        self.data = None
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

        self.b = 0
        self.d = 0
        self.c = 0
        self.minCcn = 4
        self.minDp = 0
        self.maxDp = 0
        self.maxDpAsym = 0
        self.minDpAsym = 0

        self.bList = []
        self.dList = []
        self.cList = []

        #Graph variables
        self.originalGraph = None
        self.adjustedGraph = None
        self.dryDiaGraph = None
        self.totalGraph = None

        self.adjustedGraphList = []
        self.dryDiaGraphList = []
        self.originalGraphList = []

        os.chdir(self.tempFolder)

    def setFolder(self,folder):
        """
        Set the folder where data is stored
        :param folder: path to folder
        :return:
        """
        if folder != self.folder:
            self.folder = folder
            self.run()

    def getFileNames(self):
        """
        Get all the SMPS and CCNC file names from the directory
        :return:
        """
        self.smpsTxtFilePath = []
        self.ccncCsvFilePath = []
        for root, dirs, files in os.walk(self.folder):
            for aFile in files:
                if aFile.lower().endswith('.txt'):
                    if aFile.lower().startswith("as"):
                        self.smpsTxtFilePath.append(os.path.join(self.folder, aFile))
                if aFile.lower().endswith('.csv'):
                    if aFile.lower().startswith("ccn"):
                        self.ccncCsvFilePath.append(os.path.join(self.folder,aFile))
            break

    def mergeCSVFiles(self):
        """
        Process the smps and ccnc files
        :return:
        """
        self.smpsTxtFilePath = [str(x) for x in self.smpsTxtFilePath]
        self.ccncCsvFilePath = [str(x) for x in self.ccncCsvFilePath]
        if len(self.smpsTxtFilePath) == 1:
            self.smpsFile = self.smpsTxtFilePath[0]
            print self.smpsFile
        elif len(self.smpsTxtFilePath) < 1:
            raise FileNotFoundError()
        else:
            pass
        if len(self.ccncCsvFilePath) < 1:
            raise FileNotFoundError()
        elif len(self.ccncCsvFilePath) == 1:
            self.ccncFile = self.ccncCsvFilePath[0]
        else:
            pass

    def getRawDataFromFiles(self, csvFile, txtFile):
        """
        Preprocess the data to output the necessary csv files for further processing
        :param csvFile: the path to the csv file
        :param txtFile: the path to the txt file
        """
        try:
            self.csvContent = csvProcessing(csvFile)
            self.date = self.csvContent[0]
            self.txtContent = txtProcessing(txtFile)
            self.maxPeak = len(self.txtContent[0])
            self.startTimeEntries = self.txtContent[2]
            self.endTimeEntries = self.txtContent[2][1:]
            self.timeFrame = datetime.strptime(self.endTimeEntries[0],"%I:%M:%S") - datetime.strptime(self.startTimeEntries[0],"%I:%M:%S")
            self.timeFrame = int(self.timeFrame.total_seconds())
            additionalEndTime = datetime.strptime(self.endTimeEntries[-1],"%I:%M:%S") + timedelta(seconds = self.timeFrame)
            csvEndTime = datetime.strptime(self.csvContent[1][len(self.csvContent[1])-1][0],"%I:%M:%S")
            timeDiff = additionalEndTime - csvEndTime
            additionalEndTime = additionalEndTime.time()
            timeDiff = int(timeDiff.total_seconds())
            additionalEndTime = str(additionalEndTime)
            self.endTimeEntries.append(additionalEndTime)
        except Exception:
            raise FileProcessingError()

        # Add data to the original csv data to match the number of scan with smps for further alignment
        # the last currPeak of the smps/ccnc scan should never be taken into account

        addLine = self.csvContent[1][len(self.csvContent[1])-1]
        for i in range(timeDiff):
            csvEndTime = csvEndTime + timedelta(seconds=1)
            addLine[0] = str(csvEndTime.time())
            self.csvContent[1].append(addLine)

    def getDNlog(self):
        """
        Get the DNlog data from the txt file and saves it to a csv file for easier further processing
        :param txtFile: the path to the txt file
        :return: save the date to a csv file
        """
        try:
            # Start producing the data file
            for k in range(3, len(self.txtContent)):
                if re.search('[a-zA-Z]', self.txtContent[k][0]):
                    dNlogPos = k
                    break

            startTime = ['dp'] + self.startTimeEntries
            endTime = ['dp'] + self.endTimeEntries
            dNlogList = [startTime, endTime]
            for k in range(3, dNlogPos):
                dNlogList.append(self.txtContent[k][:-1])
            self.dNlog = dNlogList
        except Exception:
            raise DNlogDataError()

    def getSMPSAndCCNC(self):
        """
        get the SMPS data from SMPS file and CCN data from CCN file
        and then combine them into a single csv file for easier further processing
        :param csvFile: the path to the csv file
        :param txtFile: the path to the txt file
        :return: print out the csv file
        """
        try:
            # Get the data from the SMPS file
            csvContent = self.csvContent[1]
            #copy the variable for local usage
            txtContent = self.txtContent
            startTime = self.startTimeEntries
            width = self.maxPeak

            # get the count position
            for k in range(3, len(self.txtContent)):
                if re.search('[a-zA-Z]', self.txtContent[k][0]):
                    for j in range(k, len(self.txtContent)):
                        if not re.search('[a-zA-Z]', self.txtContent[j][0]):
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
            sizeList = [0.625] + [0.875]
            size = 1.25
            binPos = 25
            for i in range(0, width - 1):
                sizeList.append(size)
                size += 0.5
            # get the position of the time and the average size for each scan
            for i in range(count):
                aCountList = [i + 1]
                aSizeList = [i + 1]
                for j in range(0, len(startTime)):
                    timeStamp = startTime[j]
                    for k in range(0, len(csvContent)):
                        if csvContent[k][0] == timeStamp:
                            break
                    sizeSum = 0
                    countSum = 0
                    aCountList.append(csvContent[k + i][-3])
                    for m in range(0, 20):
                        sizeSum += sizeList[m] * float(csvContent[k + i][binPos + m])
                        size += 0.5
                        countSum += float(csvContent[k + i][binPos + m])
                        # get the average size for each scan
                    if countSum == 0:
                        aSizeList.append(0)
                    else:
                        aSizeList.append(sizeSum / countSum)
                ccnList.append(aCountList)
                aveSizeList.append(aSizeList)

            # Combine the SMPS and CCNC data into the same file
            # and in a continuous form for easier processing

            # aLine has scan Time, dp, real time, SMPS count ccnc count, and ave Size
            aLine = []
            timeStamp = startTime[0]
            timeStamp = datetime.strptime(timeStamp, "%I:%M:%S")

            title = ["scan time"] + ["real time"] + ["dp"] + ["SMPS Count"] + ["CCNC Count"] + ["Ave Size"]
            smpsCcnList = []
            for i in range(0, width):
                for j in range(len(ccnList)):
                    aLine = [smpsList[j][1]] + [timeStamp.time()] + [float(smpsList[j][0])] + [float(smpsList[j][i + 2])] + [
                        float(ccnList[j][i + 1])] + [float(aveSizeList[j][i + 1])]
                    timeStamp = timeStamp + timedelta(seconds=1)
                    smpsCcnList.append(aLine)
            self.data = pandas.DataFrame(smpsCcnList, columns=title)
        except Exception:
            raise DataError()

    def matchSMPSCCNCData(self):
        """
            Calculate the stretch factor of the ccnc measurement
            :return: a number indicating the stretch factor
            """
        try:
            peakCountCCNCList = []
            peakCountSMPSList = []
            startTime = 0
            endTime = self.timeFrame
            newData = [self.data.columns.values.tolist()]



            # need to do something to calculate the correct min Dist
            minDist = 16




            shiftFactor = None
            aPeakDataForCountPeak = numpy.asarray(self.data[startTime:endTime]["SMPS Count"])
            indexes = peakutils.indexes(aPeakDataForCountPeak, thres=0.6, min_dist=minDist)
            # assuming that count currPeak of smps is always correct, the first currPeak is in the right position
            minPosSMPS = indexes[0] + getMinIndex(aPeakDataForCountPeak[indexes[0]:indexes[1]])

            # count currPeak ccnc
            aPeakDataForCountPeak = numpy.asarray(self.data[startTime:endTime]["CCNC Count"])
            indexes = peakutils.indexes(aPeakDataForCountPeak, thres=0.6, min_dist=minDist)
            minPosCCNC = indexes[0] + getMinIndex(aPeakDataForCountPeak[indexes[0]:])
            shiftFactor = minPosCCNC - minPosSMPS

            while True:
                # count currPeak smps
                aPeakDataForCountPeak = numpy.asarray(self.data[startTime:endTime]["SMPS Count"])
                indexes = peakutils.indexes(aPeakDataForCountPeak, thres=0.6, min_dist=minDist)
                # assuming that count currPeak of smps is always correct, the first currPeak is in the right position
                minPos = indexes[0] + getMinIndex(aPeakDataForCountPeak[indexes[0]:indexes[-1]])
                indexes = [x + startTime for x in indexes]
                peakCountSMPSList.append(minPos + startTime)
                # count currPeak ccnc
                aPeakDataForCountPeak = numpy.asarray(
                    self.data[startTime + shiftFactor:endTime + shiftFactor]["CCNC Count"])
                indexes = peakutils.indexes(aPeakDataForCountPeak, thres=0.6, min_dist=minDist)
                minPos = indexes[0] + getMinIndex(aPeakDataForCountPeak[indexes[0]:indexes[-1]])
                indexes = [x + startTime for x in indexes]
                peakCountCCNCList.append(minPos + startTime)
                for i in range(self.timeFrame):
                    scanTime = i + 1
                    realTime = self.data.iat[startTime + i, 1]
                    dp = self.data.iat[startTime + i, 2]
                    smpsCount = self.data.iat[startTime + i, 3]
                    ccncCount = self.data.iat[startTime + i + shiftFactor, 4]
                    aveSize = self.data.iat[startTime + i + shiftFactor, 5]
                    newData.append([scanTime, realTime, dp, smpsCount, ccncCount, aveSize])
                # loop end condition
                # use end time to check to exclude the last currPeak
                startTime += self.timeFrame
                endTime += self.timeFrame
                if endTime >= len(self.data):
                    break
            headers = newData.pop(0)
            self.data = pandas.DataFrame(newData, columns=headers)
            self.data = dateConvert(self.data)
        except:
            raise DataMatchingException()

    def finalizeData(self):
        """
        Prepare the data in the form of the shiftCalculatorSmooth to finalize the calculation for one currPeak
        :return: a list, which contains the data as in the ShiftCalculatorSmooth
        """
        ccnCNList = []
        dNdLogDpList = []
        startPoint = self.currPeak * self.timeFrame
        endPoint = (self.currPeak +1 )* self.timeFrame
        self.diameterList = list(self.data[startPoint:endPoint]['dp'])
        ccnList = list(self.data[startPoint:endPoint]['CCNC Count'])
        cnList = list(self.data[startPoint:endPoint]['SMPS Count'])

        #modify this to get a decently good number
        cnList = [x * 0.2 for x in cnList]
        self.ccnList = ccnList
        self.cnList = cnList
        self.diameterMidpointList = []
        self.ccncnList = ccnCNList
        for i in range(len(ccnList)):
            self.ccncnList.append(ccnList[i] / cnList[i])
        self.dropSizeList = list(self.data[startPoint:endPoint]['Ave Size'])
        for i in range(2,len(self.dNlog)):
            self.diameterMidpointList.append(self.dNlog[i][0])
            dNdLogDpList.append(self.dNlog[i][self.currPeak + 1])
        self.ccnNormalizedList = normalizeList(dNdLogDpList)

    def initCorrectCharges(self):
        self.cnFixedList = self.cnList[:]
        self.ccnFixedList = self.ccnList[:]
        self.gCcnList = self.ccnFixedList[:]
        self.gCnList = self.cnFixedList[:]

    def correctCharges(self):
        ##add code here to automatically identify the asymptope
        asymp = 99999
        """
        Correct the charges
        :return:
        """
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
        coeficientList = [[-0.0003,-0.1014,0.3073,-0.3372,0.1023,-0.0105],[-2.3484,0.6044,0.48,0.0013,-0.1553,0.032],[-44.4756,79.3772,-62.89,26.4492,-5.748,0.5049]]
        # frac0List = fractionCalculation(self.diameterList,0,coeficientList[0])
        frac1List = fractionCalculation(self.diameterList,1,coeficientList[1])
        frac2List = fractionCalculation(self.diameterList,2,coeficientList[2])
        frac3List = fractionCalculation(self.diameterList,3)
        chargeList = []
        for i in self.diameterList:
            QtGui.qApp.processEvents()
            aDList = [0]
            for k in range(1,4):
                c = calCC(i * 10 ** -9, lambdaAir)
                dp = 10 ** 9 * findDp(i * 10 ** -9 / c, lambdaAir, k)
                aDList.append(dp)
            chargeList.append(aDList)

        #second part of correct charges
        self.cnFixedList = self.cnList[:]
        self.ccnFixedList = self.ccnList[:]
        maxUpperBinBound = (self.diameterList[-1] + self.diameterList[-2]) / 2
        lenDpList = len(self.diameterList)
        for i in range(lenDpList):

            QtGui.qApp.processEvents()

            n = lenDpList - i - 1
            moveDoubletCounts = frac2List[n] / (frac1List[n] + frac2List[n] + frac3List[n]) * self.cnList[n]
            moveTripletCounts = frac3List[n] / (frac1List[n]  + frac2List[n]  +frac3List[n]) * self.cnList[n]
            self.cnFixedList[n] = self.cnFixedList[n] - moveDoubletCounts - moveTripletCounts
            self.ccnFixedList[n] = self.ccnFixedList[n] - moveDoubletCounts - moveTripletCounts
            if chargeList[n][2] <= maxUpperBinBound:
                j = lenDpList - 2
                while (True):
                    upperBinBound = (self.diameterList[j] +self.diameterList[j + 1]) / 2
                    lowerBinBound = (self.diameterList[j] + self.diameterList[j - 1]) / 2
                    if upperBinBound > chargeList[n][2] >= lowerBinBound:
                        self.cnFixedList[j] = self.cnFixedList[j] + moveDoubletCounts
                        if chargeList[n][2] < asymp:
                            if self.gCcnList[j] > epsilon:
                                self.ccnFixedList[j] = self.ccnFixedList[j] + moveDoubletCounts * self.gCcnList[j] / self.gCnList[j]
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
                            self.ccnFixedList[j] = self.ccnFixedList[j] + moveTripletCounts * self.ccnList[j] / self.cnList[j]
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

    def getConstants(self):
        """
        """
        # determine minDp and minDpAsym

        # The number of the data point that must produce a continuous increment to determine a vertical
        asymList = getAsym(self.diameterList, self.ccncSigList)
        increLength = int(len(asymList) / 10)
        exitLoop = False
        for i in range(len(asymList)):
            if asymList[i] > 1:
                increCount = 0
                self.minDp = self.diameterList[i]
                for j in range(i+1, i + increLength):
                    if asymList[i] > 1:
                        increCount += 1
                if increCount > increLength / 2:
                    exitLoop = True
            if exitLoop is True:
                break

        exitLoop = False
        for j in range(i + 1, len(asymList) - increLength):
            if asymList[j] < 1:
                increCount = 0
                self.minDpAsym = self.diameterList[j]
                for k in range(j + 1, j + increLength):
                    if asymList[k] < 1:
                        increCount += 1
                if increCount > increLength / 2:
                    exitLoop = True
            if exitLoop is True:
                break

        # determine maxDp
        for i in range(j, len(self.ccnFixedList)):
            if abs(self.ccncSigList[i] - self.ccncSigList[i-1]) > 0.1:
                self.maxDpAsym = self.diameterList[i]
                break
        self.maxDp = self.maxDpAsym

        asymsList = []
        for i in range(len(self.diameterList)):
            if self.minDpAsym < self.diameterList[i] < self.maxDpAsym:
                asymsList.append(self.ccncSigList[i])
            else:
                asymsList.append(0)
        self.b = getAveNoneZero(asymsList)
        self.ccncnSimList.append(0)
        for i in range (1,len(self.diameterList)):
            if self.minDp < self.diameterList[i] < self.maxDp:
                n = self.b / (1 + (self.diameterList[i] / self.d)**self.c)
                self.ccncnSimList.append(n)
            else:
                self.ccncnSimList.append(self.ccncnSimList[i - 1])

    def optimize(self):
        try:
            xList = []
            yList = []
            for i in range(len(self.diameterList)):
                if self.minDp < self.diameterList[i] < self.maxDp:
                    xList.append(self.diameterList[i])
                    yList.append(self.ccncSigList[i])

            initGuess = [30,-10]
            xList  = numpy.asarray(xList)
            yList = numpy.asarray(yList)
            # initGuess = numpy.asarray(initGuess)
            result = opt.curve_fit(f,xList, yList,bounds = ([self.minDp, -200],[self.maxDp, -1]), method = "trf")
            self.d = result[0][0]
            self.c = result[0][1]
            self.ccncnSimList = [0]
            for i in range(1, len(self.diameterList)):
                if self.minDp < self.diameterList[i] < self.maxDp:
                    n = self.b / (1 + (self.diameterList[i] / self.d) ** self.c)
                    self.ccncnSimList.append(n)
                else:
                    self.ccncnSimList.append(self.ccncnSimList[i - 1])
        except:
            raise OptimizationError()

    def makeAdjustedGraph(self):
        data  = pandas.DataFrame({"SMPS": self.cnList, "CCNC": self.ccnList})
        graph = data.plot(kind = 'line')
        graph.set_axis_bgcolor("white")
        graph.grid(color = '0.65')
        self.adjustedGraphList.append(plt.gcf())

    def makeDryDiamterGraph(self):
        plt.figure()
        plt.plot(self.diameterList,self.ccncnList, 'ro')
        plt.plot(self.diameterList, self.ccncSigList, 'bo')
        plt.plot(self.diameterMidpointList, self.ccnNormalizedList)
        plt.plot(self.diameterList, self.ccncnSimList)
        plt.gca().axes.set_ylim([-0.2, 1.2])
        plt.gca().set_axis_bgcolor("white")
        plt.grid(color = '0.65')
        self.dryDiaGraphList.append(plt.gcf())
        # graph.set_axis_bgcolor("white")
        # graph.grid(color='0.65')

    def makeGraphs(self):
        """
        Handle all of the graph making of the program
        :return:
        """
        # plt.ioff()
        self.makeAdjustedGraph()
        self.makeDryDiamterGraph()


    def preparationProcedure(self):
        try:
            self.makeProgress(maxValue = 7)
            self.makeProgress("Searching folder for CCNC and SMPS files...")
            self.getFileNames()
            self.makeProgress("Merging SMPS and CCNC files...")
            self.mergeCSVFiles()
            self.makeProgress("Processing raw data from files...")
            self.getRawDataFromFiles(self.ccncFile, self.smpsFile)
            self.makeProgress("Acquiring DnLog data...")
            self.getDNlog()
            self.makeProgress("Acquiring SMPS and CCNC data...")
            self.getSMPSAndCCNC()
            self.makeProgress("Transforming the CCNC data to match SMPS data....")
            self.matchSMPSCCNCData()
            self.makeProgress(complete = True)
        except FileNotFoundError:
            self.view.showError("Can't file SMPS or CCNC files in folder: " + str(self.folder))
            raise DataPreparationError()
        except FileProcessingError:
            self.view.showError("Can't process the SMPS or CCNC files!")
            raise DataPreparationError()
        except DNlogDataError:
            self.view.showError("Can't process DNlog data from the SMPS file!")
            raise DataPreparationError()
        except DataError:
            self.view.showError("Can't process data from SMPS or CCNC files!")
            raise DataPreparationError()
        except DataMatchingException:
            self.view.showError("Can't match SMPS and CCNC data!")
            raise DataPreparationError()

    def ProcessingProcedure(self):
        try:
            self.ccncSigList = []
            self.finalizeData()
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
            self.optimize()
            self.makeProgress()
            # Acquire data produced
            self.bList.append(self.b)
            self.dList.append(self.d)
            self.cList.append(self.c)
            self.ccnNormalizedFullList.append(self.ccnNormalizedList)
            self.ccncnFullList.append(self.ccncnList)
            self.ccncSigFullList.append(self.ccncSigList)
            self.ccncnSimFullList.append(self.ccncnSimList)
            self.makeGraphs()
            self.makeProgress()
        except OptimizationError:
            # self.view.showError("The data is not optimizable. No optimal soluation found!")
            print "There is an error on peak"


    def getPeakData(self):
        # startPoint = self.currPeak * self.timeFrame
        # endPoint = (self.currPeak + 1) * self.timeFrame
        # self.diameterList = list(self.data[startPoint:endPoint]['dp'])
        # self.ccnList = list(self.data[startPoint:endPoint]['CCNC Count'])
        # self.cnList = list(self.data[startPoint:endPoint]['SMPS Count'])
        # self.cnList = [x * 0.2 for x in self.cnList]
        # self.ccncnList = self.ccncnFullList[self.currPeak]
        # self.dropSizeList = list(self.data[startPoint:endPoint]['Ave Size'])
        # self.ccnNormalizedList = self.ccnNormalizedFullList[self.currPeak]
        self.b = self.bList[self.currPeak]
        self.d = self.dList[self.currPeak]
        self.c = self.cList[self.currPeak]
        self.adjustedGraph = self.adjustedGraphList[self.currPeak]
        self.dryDiaGraph = self.dryDiaGraphList[self.currPeak]


    def run(self):
        """
        The main running procedure of the program
        :return:
        """
        if not self.mode:
            if not self.view:
                view = MainWindow(self)
                self.setView(view)
                view.run()
            else:
                try:
                    self.preparationProcedure()
                except:
                    pass
                else:
                    # -1 because we kind of making up the last peak for alignment
                    # 5 for 5 step for each peak procesisng
                    self.makeProgress( "Preparing data for peak processing...", maxValue = self.maxPeak * 12)
                    for i in range(0, self.maxPeak - 1):
                        self.currPeak = i
                        self.makeProgress("Processing peak " + str(i + 1), value = 0)
                        # time.sleep(0.1)
                        self.ProcessingProcedure()
                    # reset currPeak
                    self.currPeak = -1
                    self.makeProgress(complete = True)
                    self.updatePeak(0)

        else:
            self.getRawDataFromFiles(csvFilePath, txtFilePath)
            self.getDNlog()
            self.getSMPSAndCCNC()
            self.matchSMPSCCNCData()
            self.finalizeData()
            removeSmallCcn(self.ccnList, self.minCcn)
            self.initCorrectCharges()
            # Since the program requires to press the button a few times
            # We are going to simulate that action by looping
            for i in range(5):
                self.correctCharges()
            for i in range(len(self.ccnFixedList)):
                self.ccncSigList.append(self.ccnFixedList[i] / self.cnFixedList[i])
            self.getConstants()
            self.optimize()
            self.makeGraphs()

    def updatePeak(self,peak):
        if peak != self.currPeak:
            self.currPeak = peak
            self.getPeakData()
            self.changeView()
            print "123"

    def changeView(self):
        self.view.updateFigures(self.adjustedGraph, self.dryDiaGraph)

    def setView(self,view):
        """
        Set the view for the controller
        :param view: the view which associate with the controller
        """
        self.view = view

    def showData(self):
        if not self.mode:
            self.view.showData()

    def makeProgress(self,message = None, maxValue = None, complete = False, value = 1):
        if not self.mode:
            self.view.makeProgress(message, maxValue, complete, value)

        else:
            print message

def removeSmallCcn(ccnList, minValue):
    for i in range(len(ccnList)):
        if ccnList[i] < minValue:
            ccnList[i] = 0

def getAveNoneZero(aList):
    sum = 0
    count = 0
    for i in aList:
        if i != 0:
            sum += i
            count += 1
    if count != 0:
        return sum / count
    else:
        return 0

def f(x, d, c):
    """
    The function for the optimization method
    :param x: the variable
    :param d: the first optimizing parameter
    :param c: the second optimizing parameter
    :return:
    """
    return 0.903 / (1+(x/d)**c)


def getAsym(xList, yList):
    asymList = []
    for i in range(0, len(xList)-1):
        if (xList[i+1] - xList[i]) != 0:
            asymList.append((yList[i+1] - yList[i])/ (xList[i+1] - xList[i]) * 100)
        else:
            asymList.append(9999)
        if 0 > asymList[i] > -0.001:
            asymList[i] = 0
    return asymList


def normalizeList(aList):
    """
    Normalize a list by divide for the largest number
    :param aList: the input list
    :return: the output list after normalized
    """
    aList = [float(x) for x in aList]
    maxValue = max(aList[10:])
    for x in range(len(aList)):
        aList[x] /= maxValue
    return aList

def csvProcessing(filePath):
    """
    read in the csv file with the name filePath
    :param filePath: the path to the csv file
    :return: a list, which represents the csv
    """
    with open(filePath,'r') as csvFile:
        reader = csv.reader(csvFile, delimiter=' ')
        line = 0
        # convert to list for easier processing
        csvContent = []
        for row in reader:
            csvContent.append(row)

        # get the date
        date = csvContent[1][1]
        date =date.replace(",", "")

        # get the time
        time = csvContent[2][0]
        time = time.replace("Time,", "")

        # clean the data for easier processing
        csvContent = filter(None, csvContent)

        # remove the unnecessary data
        csvContent = csvContent[3:]

        # process the csv
        for i in range(0,len(csvContent)):
            # remove empty string in each row
            csvContent[i] = filter(None,csvContent[i])
            # remove the "," character in each item
            for j in range(0, len(csvContent[i])):
                csvContent[i][j] = csvContent[i][j].replace(",", "")
        return date, csvContent

def txtProcessing (filePath):
    """
    read in the txt file and process the txt with the path filePath
    :param filePath: the path to the file
    :return: a list, which represents the txt
    """
    with open(filePath, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        line = 0
        # convert to list for easier processing
        txtContent = []
        for row in reader:
            txtContent.append(row)

        # remove the unnecessary head
        txtContent = txtContent[13:]

        # remove the unnecessary title
        for i in range(0,3):
            txtContent[i] = txtContent[i][1:]
        txtContent = txtContent[0:3] + txtContent[4:]

        # clean the txt file
        for i in range(0, len(txtContent)):
            txtContent[i] = filter(None, txtContent[i])
        return txtContent

def getMinIndex(data):
    """
    get the position of the smallest value of a list
    :param data: the data to process
    :return: the index of the smallest value
    """
    minValue = data[0]
    minPos = 0
    for i in range(len(data)):
        if data[i] < minValue:
            minPos = i
            minValue = data[i]
    return minPos

def dateConvert(df):
    dt = df.index
    df['real time'] = dt
    df.reset_index(drop=True)
    return df

def makeMinGraph(smpsList, ccncList):
    """
 m
    :param smpsList:
    :param ccncList:
    :return:
    """

def printList(aList):
    """
    print the list in a nice format
    :param aList: the list to print out
    :return:
    """
    for row in aList:
        print(row)
    print()

def fractionCalculation(aList, chargeNumber, coeffList = None):
    """
    :param aList: the input list
    :param coeffList: the list of coefficient
    :return: the list of fraction after calculation
    """
    newList = []
    epsilon = 0.0000001
    e = scipy.constants.e
    e0 = scipy.constants.epsilon_0
    k = scipy.constants.k
    t = scipy.constants.zero_Celsius + 25
    z = 0.875
    p = 1013
    nair = 0.000001458 * t**1.5 / (t + 110.4)
    lambdaAir = 2*nair/100/p/(8*28.84/pi/8.314/t)**0.5*1000**0.5

    if (chargeNumber <= 2):
        for i in range(len(aList)):
            l = log10(aList[i])
            sum = 0
            for j in range(len(coeffList)):
                sum += coeffList[j] * l**j
            newList.append(10 ** sum)
    else:
        for i in (aList):
            #divide the equation for easier processing
            f1 = e / sqrt(4 * pi ** 2 *e0 *i * 10 **(-9) *k * t)
            f2 = -((chargeNumber - (2 * pi *e0 *i * 10**(-9) *k *t /e ** 2)*log(z)) ** 2)
            f3 = (4 * pi *e0 *i * 10 **(-9) *k *t /e ** 2)
            value = f1 * exp(f2 / f3)
            newList.append(value)

    return newList

def calCC(dp, lambdaAir):
    return 1 + 2 * lambdaAir / dp * (1.257 + 0.4 * exp(-1.1 * dp / 2 / lambdaAir))

def findDp(dp, lambdaAir, n):
    dpOld = dp * 1 * n
    for i in range(1000):
        c = calCC(dpOld, lambdaAir)
        dpNew = dp * c * n
        if abs(dpNew - dpOld) / dpOld < 0.000001:
            break
        else:
            dpOld = dpNew
    return dpNew

def aveList(aList):
    return sum(aList) / len(aList)

def movingAve(aList, n):
    """
    Calculate a new list based on the moving average of n elements
    :param aList: the list to calculate
    :param n: the number of points to calculate for a average
    :return: a new list, which is the result of the moving average
    """
    newList = []
    for i in range(len(aList) - n + 1):
        newElement = aveList(aList[i:i + n])
        newList.append(newElement)
    return newList

def cleanseZero(aList):
    """
    Change all 0 in the list to a very small number
    :param aList: the input list with 0s
    :return: a new list with all 0s replaced
    """
    epsilon = 0.000001
    for i in range(len(aList)):
        if aList[i] == 0:
            aList[i] = epsilon
    return aList

def main():
    controller = Controller(False)
    controller.run()

if __name__ == '__main__':
    main()


