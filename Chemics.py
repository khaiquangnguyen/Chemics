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



csvFilePath = "E:\Updated program\Demonstration_Files\CCN data 100203092813.csv"
txtFilePath = "E:\Updated program\Demonstration_Files\AS_Calibration_SMPS.txt"

class Controller:
    def __init__(self,mode = True):
        """"
        The controller class for the program. ALso handle the task of the model
        """
        # True if testing, false if not
        self.ccncnSimList = []
        self.mode = mode
        # The view associated with the controller
        self.view = None
        # The folder which contains the data files
        self.folder = None
        # The list of paths to the txt files
        self.smpsTxtFilePath = []
        # The list of paths to the cccn files
        self.ccncCsvFilePath = []
        # The current analysing peak
        self.currPeak = 0
        # The total number of peaks
        self.numPeak = 0
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
        # the date of the experiment
        self.date = None
        # the time frame of a a measurement
        self.timeFrame = 135
        #the content of the csv file
        self.csvContent = None
        #the content of the txt file
        self.txtContent = None
        #The list of all diameter of measurement
        self.diameterList = None
        #The DNlog data
        self.dNlog = None
        #The SMPS and CCNC data
        self.data = None
        #CCN List
        self.ccnList = None
        #CN List
        self.cnList = None
        self.ccnFixedList = None
        self.cnFixedList = None
        self.gCnList = None
        self.gCcnList = None
        self.ccncnSigList = []
        self.ccnNormalizedList = []
        self.diameterMidpointList = []
        self.asyms = []
        self.b = 0
        self.d = 33.3667106871589
        self.c  = -18.3807247490387
        self.minCcn = 4

        # os.chdir(self.tempFolder)

    def setFolder(self,folder):
        """
        Set the folder where data is stored
        :param folder: path to folder
        :return:
        """
        if folder != self.folder:
            self.folder = folder
            self.getFileNames()

    def getFileNames(self):
        """
        Get all the SMPS and CCNC file names from the directory
        :return:
        """
        for root, dirs, files in os.walk(self.folder):
            for aFile in files:
                if aFile.lower().endswith('.txt'):
                    if aFile.lower().startswith("as"):
                        self.smpsTxtFilePath.append(os.path.join(self.folder, aFile))
                if aFile.lower().endswith('.csv'):
                    if aFile.lower().startswith("ccn"):
                        self.ccncCsvFilePath.append(os.path.join(self.folder,aFile))
        self.mergeCSVFiles()

    def mergeCSVFiles(self):
        """
        Process the smps and ccnc files
        :return:
        """
        self.smpsTxtFilePath = [str(x) for x in self.smpsTxtFilePath]
        self.ccncCsvFilePath = [str(x) for x in self.ccncCsvFilePath]
        if len(self.smpsTxtFilePath) !=1:
            pass
        if len(self.ccncCsvFilePath) != 1:
            pass
        self.smpsFile = self.smpsTxtFilePath[0]
        self.ccncFile = self.ccncCsvFilePath[0]

    def getRawDataFromFiles(self, csvFile, txtFile):
        """
        Preprocess the data to output the necessary csv files for further processing
        :param csvFile: the path to the csv file
        :param txtFile: the path to the txt file
        """
        self.csvContent = csvProcessing(csvFile)
        self.date = self.csvContent[0]
        self.makeProgress("Finish reading csv. Reading txt ...")
        self.txtContent = txtProcessing(txtFile)
        self.makeProgress("Finish reading txt. Getting necessary information...")
        self.numPeak = len(self.txtContent[0])
        self.startTimeEntries = self.txtContent[2]
        self.endTimeEntries = self.txtContent[2][1:]


        # get the number of count
        loopCount = 0
        additionalEndTime = datetime.strptime(self.endTimeEntries[-1],"%I:%M:%S") + timedelta(seconds = self.timeFrame)
        csvEndTime = datetime.strptime(self.csvContent[1][len(self.csvContent[1])-1][0],"%I:%M:%S")
        timeDiff = additionalEndTime - csvEndTime
        additionalEndTime = additionalEndTime.time()
        timeDiff = int(timeDiff.total_seconds())
        additionalEndTime = str(additionalEndTime)
        self.endTimeEntries.append(additionalEndTime)

        #add data to the original csv data to match the number of scan with smps for further alignment
        #the last peak of the smps/ccnc scan should never be taken into account

        addLine = self.csvContent[1][len(self.csvContent[1])-1]
        for i in range(timeDiff):
            csvEndTime = csvEndTime + timedelta(seconds=1)
            addLine[0] = str(csvEndTime.time())
            self.csvContent[1].append(addLine)
        self.makeProgress("Finish getting necessary information. Begin getting data...")
        self.makeProgress("Finish getting SMPS and CCNC dp data....")
        self.makeProgress("Finish processing raw data. Working on combined data...")

    def getDNlog(self):
        """
        get the DNlog data from the txt file and saves it to a csv file for easier further processing
        :param txtFile: the path to the txt file
        :return: save the date to a csv file
        """

        # start producing the data file
        # get the dNlog position
        for k in range(3, len(self.txtContent)):
            if re.search('[a-zA-Z]', self.txtContent[k][0]):
                dNlogPos = k
                break

        startTime = ['dp'] + self.startTimeEntries
        endTime = ['dp'] + self.endTimeEntries
        dNlogList = [startTime, endTime]
        for k in range(3, dNlogPos):
            dNlogList.append(self.txtContent[k][:-1])

        # Write the data to csv file
        with open('dNlog.csv', "wb") as f:
            writer = csv.writer(f)
            writer.writerows(dNlogList)
        self.dNlog = dNlogList
        self.makeProgress("Finish getting DNlog data. Beging getting Dp data...")

    def getSMPSAndCCNC(self):
        """
        get the SMPS data from SMPS file and CCN data from CCN file
        and then combine them into a single csv file for easier further processing
        :param csvFile: the path to the csv file
        :param txtFile: the path to the txt file
        :return: print out the csv file
        """
        # Get the data from the SMPS file
        csvContent = self.csvContent[1]
        #copy the variable for local usage
        txtContent = self.txtContent
        startTime = self.startTimeEntries
        width = self.numPeak

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
        smpsCcnList = [title]
        for i in range(0, width):
            for j in range(len(ccnList)):
                aLine = [smpsList[j][1]] + [timeStamp.time()] + [smpsList[j][0]] + [smpsList[j][i + 2]] + [
                    ccnList[j][i + 1]] + [aveSizeList[j][i + 1]]
                timeStamp = timeStamp + timedelta(seconds=1)
                smpsCcnList.append(aLine)

        with open('smps_ccn.csv', "wb") as f:
            writer = csv.writer(f)
            writer.writerows(smpsCcnList)

    def matchSMPSCCNCData(self):
        """
        Calculate the stretch factor of the ccnc measurement
        :param data: the data
        :return: a number indicating the stretch factor
        """
        data = pandas.read_csv("smps_ccn.csv")
        peakCountCCNCList = []
        peakCountSMPSList = []
        startTime = 0
        endTime = self.timeFrame
        newData = [data.columns.values.tolist()]
        # need to do something to calculate the correct min Dist
        minDist = 16
        shiftFactor = None
        aPeakDataForCountPeak = numpy.asarray(data[startTime:endTime]["SMPS Count"])
        indexes = peakutils.indexes(aPeakDataForCountPeak, thres=0.6, min_dist=minDist)
        # assuming that count peak of smps is always correct, the first peak is in the right position
        minPosSMPS = indexes[0] + getMinIndex(aPeakDataForCountPeak[indexes[0]:indexes[1]])

        # count peak ccnc
        aPeakDataForCountPeak = numpy.asarray(data[startTime:endTime]["CCNC Count"])
        indexes = peakutils.indexes(aPeakDataForCountPeak, thres=0.6, min_dist=minDist)
        minPosCCNC = indexes[0] + getMinIndex(aPeakDataForCountPeak[indexes[0]:])



        # shiftFactor = minPosCCNC - minPosSMPS


        shiftFactor = 8

        while True:
            # count peak smps
            aPeakDataForCountPeak = numpy.asarray(data[startTime:endTime]["SMPS Count"])
            indexes = peakutils.indexes(aPeakDataForCountPeak, thres=0.6, min_dist=minDist)
            # assuming that count peak of smps is always correct, the first peak is in the right position
            minPos = indexes[0] + getMinIndex(aPeakDataForCountPeak[indexes[0]:indexes[-1]])
            indexes = [x + startTime for x in indexes]
            peakCountSMPSList.append(minPos + startTime)
            # count peak ccnc
            aPeakDataForCountPeak = numpy.asarray(
                data[startTime + shiftFactor:endTime + shiftFactor]["CCNC Count"])
            indexes = peakutils.indexes(aPeakDataForCountPeak, thres=0.6, min_dist=minDist)
            minPos = indexes[0] + getMinIndex(aPeakDataForCountPeak[indexes[0]:indexes[-1]])
            indexes = [x + startTime for x in indexes]
            peakCountCCNCList.append(minPos + startTime)
            for i in range(self.timeFrame):
                scanTime = i + 1
                realTime = data.iloc[startTime + i, 1]
                dp = data.iloc[startTime + i, 2]
                smpsCount = data.iloc[startTime + i, 3]
                ccncCount = data.iloc[startTime + i + shiftFactor, 4]
                aveSize = data.iloc[startTime + i + shiftFactor, 5]
                newData.append([scanTime, realTime, dp, smpsCount, ccncCount, aveSize])
            # loop end condition
            # use end time to check to exclude the last peak
            startTime += self.timeFrame
            endTime += self.timeFrame
            if endTime >= len(data):
                break

        # print len(peakCountSMPSList)
        # print len(peakCountCCNCList)
        # print peakCountSMPSList
        # print peakCountCCNCList
        self.data = newData
        headers = self.data.pop(0)
        self.data = pandas.DataFrame(self.data, columns=headers)
        self.data = dateConvert(self.data)

    def makeOriginalGraph(self):
        """
         Visualize the data
        """
        self.makeProgress("Finish aligning SMPS and CCNC data...")
        dataCopied = self.data
        data = pandas.melt(self.data, id_vars="real time", value_vars=["CCNC Count", "SMPS Count"])
        fullGraph = ggplot(data, aes('real time', 'value', color='variable')) + theme_bw() + geom_step()
        fullGraph.save(fullGraphFileName, width=90, height=3.9)
        self.makeProgress("Finish producing full scale makeOriginalGraph...")

        # do not output the png if in testing mode
        if self.mode:
            return
        for i in range(self.numPeak - 1):
            startPos = i * self.timeFrame
            endPos = (i + 1) * self.timeFrame
            peakData = pandas.melt(dataCopied[startPos:endPos], id_vars="real time", value_vars=["CCNC Count", "SMPS Count"])
            fullGraph = ggplot(peakData,
                               aes('real time', 'value', color='variable')) + theme_bw() + geom_step()
            fileName = singlePeakFileName + str(i) + ".png"
            fullGraph.save(fileName, width=10, height=6)
        self.makeProgress("Finish producing individual peak makeOriginalGraph...")
        self.showData()

    def finalizeData(self):
        """
        Prepare the data in the form of the shiftCalculatorSmooth to finalize the calculation for one peak
        :return: a list, which contains the data as in the ShiftCalculatorSmooth
        """
        self.currPeak = 12
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
        self.ccncnList = ccnCNList
        for i in range(len(ccnList)):
            self.ccncnList.append(ccnList[i] / cnList[i])
        dropSizeList = list(self.data[startPoint:endPoint]['Ave Size'])
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
        epsilon = 0.0000001
        e = scipy.constants.e
        e0 = scipy.constants.epsilon_0
        k = scipy.constants.k
        t = scipy.constants.zero_Celsius + 25
        z = 0.875
        p = 1013
        nair = 0.000001458 * t ** 1.5 / (t + 110.4)
        lambdaAir = 2 * nair / 100 / p / (8 * 28.84 / pi / 8.314 / t) ** 0.5 * 1000 ** 0.5
        coeficientList = [[-0.0003,-0.1014,0.3073,-0.3372,0.1023,-0.0105],[-2.3484,0.6044,0.48,0.0013,-0.1553,0.032],[-44.4756,79.3772,-62.89,26.4492,-5.748,0.5049]]
        frac0List = fractionCalculation(self.diameterList,0,coeficientList[0])
        frac1List = fractionCalculation(self.diameterList,1,coeficientList[1])
        frac2List = fractionCalculation(self.diameterList,2,coeficientList[2])
        frac3List = fractionCalculation(self.diameterList,3)
        chargeList = []
        for i in self.diameterList:
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
                    if chargeList[n][2] < upperBinBound and chargeList[n][2] >= lowerBinBound:
                        self.cnFixedList[j] = self.cnFixedList[j] + moveDoubletCounts
                        if chargeList[n][2] < asymp:
                            if self.gCcnList[j] > 0.0000000001:
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
                    if chargeList[n][3] < upperBinBound and chargeList[n][3] >= lowerBinBound:
                        self.cnFixedList[j] = self.cnFixedList[j] + moveTripletCounts
                        if chargeList[n][3] < asymp:
                            self.ccnFixedList[j] = self.ccnFixedList[j] + moveTripletCounts * self.ccnList[j] / self.cnList[j]
                        else:
                            self.ccnFixedList[j] = self.ccnFixedList[j] + moveTripletCounts
                        break
                    j -= 1

        for i in range(len(self.ccnFixedList)):
            if self.ccnFixedList[i] / self.cnFixedList[i] < -0.01:
                self.ccnFixedList[i] = 0

        self.gCcnList = self.ccnFixedList[:]
        self.gCnList = self.cnFixedList[:]
        # print self.cnList
        # print cnFixedList
        # print self.ccnList
        # print ccnFixedList

    def getConstants(self):
        """
        Get the necessary constants for the sigmodal fit
        :return:
        """
        # determine minDp and minDpAsym
        increCount = 0
        #the percentage of the data that must produce a continuous increment to determine a vertical
        increLength = 0.06 * len(self.ccncnSigList)
        print len(self.ccnFixedList)
        for i in range(1,len(self.ccnFixedList)):
            if self.ccncnSigList[i] >= self.ccncnSigList[i-1]:
                if increCount == 0:
                    self.minDp = self.diameterList[i - 1]
                increCount += 1
                if increCount >= increLength:
                    for j in range(i+1, len(self.ccncnSigList)):
                        if self.ccncnSigList[j] <= self.ccncnSigList[j-1]:
                            self.minDpAsym = self.diameterList[j]
                            break
                    break

            else:
                increCount = 0
                self.minDp = 0
                self.minDpAsym = 0

        # determine maxDp
        maxCcnCn = max(self.ccncnSigList)
        limit = 0.1 * maxCcnCn
        for i in range(1, len(self.ccnFixedList)):
            if abs(self.ccncnSigList[i] - self.ccncnSigList[i-1]) > limit and i > len(self.ccncnSigList) * 2 / 3:
                self.maxDpAsym = self.diameterList[i]
                break
        self.maxDp = self.maxDpAsym

        for i in range(len(self.diameterList)):
            if self.diameterList[i] > self.minDpAsym and self.diameterList[i] < self.maxDpAsym:
                self.asyms.append(self.ccncnSigList[i])
            else:
                self.asyms.append(0)

        self.b = getAveNoneZero(self.asyms)
        self.ccncnSimList.append(0)
        for i in range (1,len(self.diameterList)):
            if self.diameterList[i] > self.minDp and self.diameterList[i] < self.maxDp:
                n = self.b / (1 + (self.diameterList[i] / self.d)**self.c)
                self.ccncnSimList.append(n)
            else:
                self.ccncnSimList.append(self.ccncnSimList[i - 1])

        print self.b

    def makeAdjustedGraph(self):
        data  = pandas.DataFrame({"SMPS": self.cnList, "CCNC": self.ccnList})
        graph = data.plot(kind = 'line')
        graph.set_axis_bgcolor("white")
        graph.grid(color = '0.65')
        plt.show()

    def makeDryDiamterGraph(self):
        # data = pandas.DataFrame({"scan time": scanList, "SMPS": self.cnList, "CCNC": self.ccnList})
        # data = pandas.melt(data, id_vars="scan time", value_vars=["SMPS", "CCNC"])
        # fullGraph = ggplot(data, aes('scan time', 'value', color='variable')) + theme_bw() + geom_step()
        # print fullGraph
        dataRaw = pandas.DataFrame({"dp": self.diameterList, "ccn/cn": self.ccncnList, "ccn/cn Fixed": self.ccncnSigList})
        # graph = dataRaw.plot(kind ='scatter', x = 'dp', y = 'ccn/cn')
        # graph = dataRaw.plot(kind='scatter', x='dp', y='ccn/cn Fixed')

        plt.plot(self.diameterList,self.ccncnList, marker = 'o')
        plt.plot(self.diameterList, self.ccncnSigList, marker='o')
        plt.plot(self.diameterMidpointList, self.ccnNormalizedList)

        # graph.set_axis_bgcolor("white")
        # graph.grid(color='0.65')
        plt.show()

    def makeGraphs(self):
        """
        Handle all of the graph making of the program
        :return:
        """
        # self.makeOriginalGraph()
        # self.makeAdjustedGraph()
        self.makeDryDiamterGraph()

    def run(self):
        """
        The main running procedure of the program
        :return:
        """
        if not self.mode:
            pass
            # if self.view:
            #     self.getRawDataFromFiles(self.ccncCsvFilePath, self.smpsTxtFilePath)
            #     self.getDNlog()
            #     self.getSMPSAndCCNC()
            #     self.matchSMPSCCNCData()
            # else:
            #     view = MainWindow(self)
            #     self.setView(view)
            #     view.run()
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
            for i in range(10):
                self.correctCharges()
            for i in range(len(self.ccnFixedList)):
                self.ccncnSigList.append(self.ccnFixedList[i] / self.cnFixedList[i])
            self.getConstants()
            # self.makeGraphs()


    def setView(self,view):
        """
        Set the view for the controller
        :param view: the view which associate with the controller
        """
        self.view = view

    def showData(self):
        if not self.mode:
            self.view.showData()

    def makeProgress(self,message = None):
        if not self.mode:
            self.view.showProgress()
            self.view.makeProgress(message)

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
    return sum / count

def normalizeList(aList):
    """
    Normalize a list by divide for the largest number
    :param aList: the input list
    :return: the output list after normalized
    """
    aList = [float(x) for x in aList]
    maxValue = max(aList)
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
    controller = Controller(True)
    controller.run()

if __name__ == '__main__':
    main()
