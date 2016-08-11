import csv
import re
import pandas
from datetime import *
import numpy
import peakutils
from scipy import *
from scipy import stats
import tempfile
import scipy.constants
from sys import exit
import settings
import time

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
    return 0.903 / (1 + (x / d) ** c)

def getAsym(xList, yList):
    asymList = []
    # Normalize xList
    xList = normalizeList(xList)
    # Normalize yList
    yList = normalizeList(yList)
    # Get the asymtope
    for i in range(0, len(xList) - 1):
        if (xList[i + 1] - xList[i]) != 0:
            asymList.append((yList[i + 1] - yList[i]) / (xList[i + 1] - xList[i]))
        else:
            asymList.append(9999)
        if 0 > asymList[i] > -0.001:
            asymList[i] = 0

    return asymList


def normalizeList(aList):
    """
    Normalize aParam list by divide for the largest number
    :param aList: the input list
    :return: the output list after normalized
    """
    # convert to float
    aList = [float(x) for x in aList]
    # remove invalids
    for i in range(len(aList)):
        if pandas.isnull(aList[i]):
            aList[i] = 0
    maxValue = max(aList[10:])
    if maxValue == 0:
        return aList
    for x in range(len(aList)):
        if aList[x]:
            aList[x] /= maxValue
    return aList


def csvProcessing(filePath):
    """
    process all csv files input
    :param filePath: the path to the csv file
    :return: aParam list, which represents the csv
    """
    if len(filePath) < 1:
        raise FileNotFoundError()
    elif len(filePath) == 1:
        return singleCSVFileProcessing(filePath[0])
    else:
        date = None
        csvContent = None
        for i in filePath:
            aCSV = singleCSVFileProcessing(i)
            date = aCSV[0]
            if csvContent:
                csvContent.extend(aCSV[1])
            else:
                csvContent = aCSV[1]
    return date, csvContent

def singleCSVFileProcessing(filePath):
    """
       read in a single csv file with path filePath
       :param filePath: the path to the csv file
       :return: aParam list, which represents the csv
       """
    with open(filePath, 'r') as csvFile:
        reader = csv.reader(csvFile, delimiter=' ')
        line = 0
        # convert to list for easier processing
        csvContent = []
        for row in reader:
            csvContent.append(row)

        # get the date
        date = csvContent[1][1]
        date = date.replace(",", "")

        # get the time
        time = csvContent[2][0]
        time = time.replace("Time,", "")

        # clean the data for easier processing
        csvContent = filter(None, csvContent)

        # remove the unnecessary data
        csvContent = csvContent[4:]

        # process the csv
        for i in range(0, len(csvContent)):
            # remove empty string in each row
            csvContent[i] = filter(None, csvContent[i])
            # remove the "," character in each item
            for j in range(0, len(csvContent[i])):
                csvContent[i][j] = csvContent[i][j].replace(",", "")
        return date, csvContent

def txtProcessing(filePath):
    """
    read in the txt file and process the txt with the path filePath
    :param filePath: the path to the file
    :return: aParam list, which represents the txt
    """
    with open(filePath, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        line = 0
        # convert to list for easier processing
        txtContent = []
        for row in reader:
            txtContent.append(row)

        # remove the unnecessary header
        for i in range(len(txtContent)):
            if ''.join(txtContent[i][0].split()).lower() == "starttime":
                txtContent = txtContent[i:]
                break
        txtContent[0] = txtContent[0][1:]
        txtContent = txtContent[0:1] + txtContent[2:]
        # clean the txt file
        for i in range(0, len(txtContent)):
            txtContent[i] = filter(None, txtContent[i])
        return txtContent

def getMinIndex(aList, threshold = 5):
    """
    get the position of the smallest value of aParam list
    :param data: the data to process
    :return: the index of the smallest value. -1 if the peak is not usable
    """
    firstMax = 0
    maxPos = 0
    maxPos = peakutils.indexes(aList,thres = 0.5, min_dist=len(aList) / 10)
    if len(maxPos) == 0:
        return -1
    else:
        maxPos = maxPos[0]
    firstMax = aList[maxPos]

    # Get the second maximum
    secondMax = 0
    secondMaxDis = 0
    secondMaxPos = 0
    for i in range(maxPos + 2, len(aList)):
            maxDis = 0
            if aList[i] < firstMax * 0.1:
                continue
            for j in range(1,i - maxPos):
                if aList[i] <= aList[i-j]:
                    break
                else:
                    maxDis += 1
            if maxDis >= secondMaxDis:
                secondMax = aList[i]
                secondMaxPos = i
                secondMaxDis = maxDis

    # Check if the two peaks are actually usable
    # Get the minimum between two peaks
    min = firstMax
    minPos = maxPos
    for i in range(maxPos, secondMaxPos):
        if aList[i] <= min:
            min = aList[i]
            minPos = i
        if aList[i]  - min < 10:
            min = aList[i]
            minPos = i


    # If the second peak is too small,then not valid
    if secondMax < firstMax * 0.1:
        return -1
    # if the minimum is too big, the not valid as well
    if min >= firstMax * 0.9:
        return -1
    # If either peak is 0
    if secondMax == 0 or firstMax == 0:
        return -1
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
    print the list in aParam nice format
    :param aList: the list to print out
    :return:
    """
    for row in aList:
        print(row)
    print()

def fractionCalculation(aList, chargeNumber, coeffList=None):
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
    nair = 0.000001458 * t ** 1.5 / (t + 110.4)
    lambdaAir = 2 * nair / 100 / p / (8 * 28.84 / pi / 8.314 / t) ** 0.5 * 1000 ** 0.5

    if (chargeNumber <= 2):
        for i in range(len(aList)):
            l = log10(aList[i])
            sum = 0
            for j in range(len(coeffList)):
                sum += coeffList[j] * l ** j
            newList.append(10 ** sum)
    else:
        for i in (aList):
            # divide the equation for easier processing
            f1 = e / sqrt(4 * pi ** 2 * e0 * i * 10 ** (-9) * k * t)
            f2 = -((chargeNumber - (2 * pi * e0 * i * 10 ** (-9) * k * t / e ** 2) * log(z)) ** 2)
            f3 = (4 * pi * e0 * i * 10 ** (-9) * k * t / e ** 2)
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
    Calculate aParam new list based on the moving average of n elements
    :param aList: the list to calculate
    :param n: the number of points to calculate for aParam average
    :return: aParam new list, which is the result of the moving average
    """
    newList = []
    for i in range(len(aList) - n + 1):
        newElement = aveList(aList[i:i + n])
        newList.append(newElement)
    return newList

def cleanseZero(aList):
    """
    Change all 0 in the list to aParam very small number
    :param aList: the input list with 0s
    :return: aParam new list with all 0s replaced
    """
    epsilon = 0.000001
    for i in range(len(aList)):
        if aList[i] == 0:
            aList[i] = epsilon
    return aList

def getCorrectNum(aList, number, bigger=True):
    """
    Get aParam number approximately around number
    Assuming the aList is sorted in descending order
    :param bigger: if bigger, return the minimum number bigger than number, otherwise
                    return the maximum number smaller than number
    :param aList: the list
    :param number: the number
    :return: the correct number
    """
    num = aList[0]
    for i in range(1, len(aList)):
        if aList[i] < number:
            if bigger:
                return (num, i - 1)
            else:
                return (aList[i], i)
        else:
            num = aList[i]
    return (aList[-1], len(aList) -1)









