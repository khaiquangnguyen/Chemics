import csv
import re
import pandas
from datetime import *
import numpy
import numpy as np
import peakutils
from scipy import *
from scipy import stats
import tempfile
import scipy.constants
from sys import exit
import settings as CONST
import time
import matplotlib.pyplot as plt
import threading
from PySide.QtCore import *
import traceback, sys


class WorkerSignals(QObject):
    started = Signal(str)
    finished = Signal()
    error = Signal(tuple)
    progress = Signal(int)


class Worker(QRunnable):
    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()
        kwargs['progress_start'] = self.signals.started
        kwargs['progress_update'] = self.signals.progress

    def run(self):
        try:
            self.fn(*self.args, **self.kwargs)
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        finally:
            self.signals.finished.emit()  # Done


def run_new_thread(fn_name, started_fn, finished_fn, progress_fn):
    worker = Worker(fn_name)
    worker.signals.started.connect(started_fn)
    worker.signals.progress.connect(progress_fn)
    worker.signals.finished.connect(finished_fn)
    return worker


def get_layout_widgets(layout):
    return (layout.itemAt(i) for i in range(layout.count()))


def smooth(a_list, method="Savitzky-Golay filter"):
    if method == "Savitzky-Golay filter":
        try:
            a_list = scipy.signal.savgol_filter(a_list, 5, 2)
        except:
            pass
    return a_list

def heavy_smooth(a_list):
    try:
        a_list = scipy.signal.savgol_filter(a_list, 15, 2)
    except:
        pass
    return a_list

def fill_zeros_to_end(a_list, length_to_fill):
    if len(a_list) < length_to_fill:
        filler_array = numpy.asarray([0] * (length_to_fill - len(a_list)))
        a_list = numpy.append(a_list, filler_array)
    return numpy.asarray(a_list)


def fill_zeros_to_begin(a_list, fill_amount):
    filler_array = numpy.asarray([0] * fill_amount)
    a_list = numpy.append(filler_array, a_list)
    return numpy.asarray(a_list)


def resolve_small_ccnc_vals(ccnc_vals):
    for i in range(len(ccnc_vals)):
        if ccnc_vals[i] < 4:
            ccnc_vals[i] = CONST.EPSILON
    return ccnc_vals


def are_floats_equal(a, b, err=0.01):
    return abs(a - b) < err


def get_ave_none_zero(a_list):
    sum = 0
    count = 0
    for i in a_list:
        if i != 0:
            sum += i
            count += 1
    if count != 0:
        return sum / count
    else:
        return 0


# def f(x, d, c):
#     """
#     The function for the optimization method
#     :param x: the variable
#     :param d: the first optimizing parameter
#     :param c: the second optimizing parameter
#     :return: a number
#     """
#     return 0.923 / (1 + (x / d) ** c)

def outliers_iqr_ver_2(value_array, index_array):
    quartile_1, quartile_3 = numpy.percentile(value_array, [25, 75])
    iqr = quartile_3 - quartile_1
    lower_bound = quartile_1 - (iqr * 1.25)
    upper_bound = quartile_3 + (iqr * 1.25)
    indexes = []
    for i in range(len(value_array)):
        if lower_bound <= value_array[i] <= upper_bound:
            indexes.append(index_array[i])
    return indexes

def get_asym_list(x_list, y_list):
    asym_list = []
    # Normalize x_list
    x_list = normalize_dndlogdp_list(x_list)
    # Normalize y_list
    y_list = normalize_dndlogdp_list(y_list)
    # Get the asymptote
    for i in range(0, len(x_list) - 1):
        if (x_list[i + 1] - x_list[i]) != 0:
            asym_list.append((y_list[i + 1] - y_list[i]) / (x_list[i + 1] - x_list[i]))
        else:
            asym_list.append(9999)
        if 0 > asym_list[i] > -0.001:
            asym_list[i] = 0
    return asym_list

def cal_percentage_more(a_list,value, min_diff):
    count = 0
    for a_val in a_list:
        if a_val > value or abs(value-a_val) < min_diff:
            count += 1
    return safe_div(count,len(a_list))

def cal_percentage_less(a_list,value, min_diff):
    count = 0
    for a_val in a_list:
        if value>a_val or abs(value-a_val) < min_diff:
            count += 1
    return safe_div(count, len(a_list))


def safe_div(x, y):
    return 0 if y == 0 else x / y

def safe_div_array(l1,l2):
    length = min(len(l1), len(l2))
    result_list = []
    for i in range(length):
        result_list.append(safe_div(l1[i],l2[i]))
    return  result_list

def create_size_list():
    # Steps to get the CCNC count. These values are hard-coded. Don't change them!
    size_list = [0.625] + [0.875]
    size = 1.25
    # Calculate the size of each bin
    for i in range(0, 18):
        size_list.append(size)
        size += 0.5
    return size_list


def normalize_dndlogdp_list(a_list):
    # the assumption here is that the first few data points are probably not quite right
    # so we only normalized based on the max after a certain point
    max_value = max(a_list[5:])
    if max_value == 0:
        return a_list
    for x in range(len(a_list)):
        if a_list[x]:
            a_list[x] /= max_value
    return numpy.asarray(a_list)


def process_csv_files(file_path):
    """
    process all csv files input
    :param file_path: the path to the csv file
    :return: aParam list, which represents the csv
    """
    if len(file_path) < 1:
        raise FileNotFoundError()
    elif len(file_path) == 1:
        return process_a_csv(file_path[0])
    else:
        date = None
        csv_content = None
        for i in file_path:
            csv = process_a_csv(i)
            date = csv[0]
            if csv_content:
                csv_content.extend(csv[1])
            else:
                csv_content = csv[1]
    return date, csv_content


def process_a_csv(file_path):
    """
       read in a single csv file with path file_path
       :param file_path: the path to the csv file
       :return: aParam list, which represents the csv
       """
    with open(file_path, 'r') as csvFile:
        reader = csv.reader(csvFile, delimiter=',')
        line = 0
        # convert to list for easier processing
        csv_content = []
        for row in reader:
            csv_content.append(row)

        # get the experiment_date
        date = csv_content[1][1]
        date = date.replace(",", "")

        # get the time
        time = csv_content[2][0]
        time = time.replace("Time,", "")

        # clean the processed_data for easier processing
        csv_content = filter(None, csv_content)

        # remove the unnecessary processed_data
        csv_content = csv_content[4:]

        # process the csv
        for i in range(0, len(csv_content)):
            # remove empty string in each row
            csv_content[i] = filter(None, csv_content[i])
            # remove the "," character in each item
            for j in range(0, len(csv_content[i])):
                csv_content[i][j] = csv_content[i][j].replace(",", "")
        return date, csv_content


def process_text_files(file_path):
    """
    read in the txt file and process the txt with the path filePath
    :param file_path: the path to the file
    :return: aParam list, which represents the txt
    """
    with open(file_path, 'r') as csv_file:
        reader = csv.reader(csv_file, delimiter='\t')
        line = 0
        # convert to list for easier processing
        txt_content = []
        for row in reader:
            txt_content.append(row)

        # remove the unnecessary header
        for i in range(len(txt_content)):
            if ''.join(txt_content[i][0].split()).lower() == "starttime":
                txt_content = txt_content[i:]
                break
        txt_content[0] = txt_content[0][1:]
        txt_content = txt_content[0:1] + txt_content[2:]
        # clean the txt file
        for i in range(0, len(txt_content)):
            txt_content[i] = filter(None, txt_content[i])
        return txt_content


def find_ref_index_smps(smps_data, up_time):
    """
    Find the reference point of smps data, based on the information that the first peak of smps data is
    between index 0 and index up_time, and the second peak is between up_time and the last index
    :param smps_data:
    :param up_time: the position that seperates the first peak and second peak
    :return: -1 if can't find min. otherwise, return the reference point
    """
    # seperate the list which contains left peak and right peak
    left_list = smps_data[0:up_time]
    right_list = smps_data[up_time:]
    # find the left peak and right peak
    left_max = numpy.argmax(left_list)
    right_max = up_time + numpy.argmax(right_list)
    # get the data between the peaks
    middle_list = smps_data[left_max:right_max]
    # then perform agressive smoothing
    middle_list = scipy.signal.savgol_filter(middle_list, 11, 2)
    # find the min between the two peaks, which is exactly what we need
    potential_mins = scipy.signal.argrelmin(middle_list, order=4)
    if len(potential_mins[0]) == 0:
        return None
    ref_index = left_max + potential_mins[0][0]
    # The reference point should be within a reasonable range of the up_time
    if abs(up_time - ref_index) > len(smps_data) / 20:
        return None
    else:
        return ref_index


def find_ref_index_ccnc(ccnc_list, potential_loc):
    # we are working with a lot of assumptions here
    # The first one is that the two ref values must be quite close to each other
    # find absolute max
    # the potential location is the sum of the smps local minimum and the base shift factor calculated from other
    # scans
    ccnc_list = smooth(ccnc_list)
    first_point = ccnc_list[potential_loc - 1]
    second_point = ccnc_list[potential_loc]
    slope = second_point - first_point
    if slope > max(ccnc_list) / 20:
        return potential_loc - 1
    else:
        return potential_loc


def convert_date(df):
    dt = df.index
    df['real time'] = dt
    df.reset_index(drop=True)
    return df


def print_list(a_list):
    """
    print the list in aParam nice format
    :param a_list: the list to print out
    :return:
    """
    for row in a_list:
        print(row)
    print()


def calculate_fraction(a_list, charge_number, coef_list=None):
    """
    :param a_list: the input list
    :param coef_list: the list of coefficient
    :return: the list of fraction after calculation
    """
    new_list = []
    epsilon = 0.0000001
    e = scipy.constants.e
    e0 = scipy.constants.epsilon_0
    k = scipy.constants.k
    t = scipy.constants.zero_Celsius + 25
    z = 0.875
    p = 1013
    nair = 0.000001458 * t ** 1.5 / (t + 110.4)
    lambda_air = 2 * nair / 100 / p / (8 * 28.84 / pi / 8.314 / t) ** 0.5 * 1000 ** 0.5

    if (charge_number <= 2):
        for i in range(len(a_list)):
            l = log10(a_list[i])
            sum = 0
            for j in range(len(coef_list)):
                sum += coef_list[j] * l ** j
            new_list.append(10 ** sum)
    else:
        for i in (a_list):
            # divide the equation for easier processing
            f1 = e / sqrt(4 * pi ** 2 * e0 * i * 10 ** (-9) * k * t)
            f2 = -((charge_number - (2 * pi * e0 * i * 10 ** (-9) * k * t / e ** 2) * log(z)) ** 2)
            f3 = (4 * pi * e0 * i * 10 ** (-9) * k * t / e ** 2)
            value = f1 * exp(f2 / f3)
            new_list.append(value)

    return new_list


def cal_cc(dp, lambda_air):
    return 1 + 2 * lambda_air / dp * (1.257 + 0.4 * exp(-1.1 * dp / 2 / lambda_air))


def find_dp(dp, lambda_air, n):
    dp_old = dp * 1 * n
    for i in range(1000):
        c = cal_cc(dp_old, lambda_air)
        dp_new = dp * c * n
        if abs(dp_new - dp_old) / dp_old < 0.000001:
            break
        else:
            dp_old = dp_new
    return dp_new


def calculate_ave_list(a_list):
    return sum(a_list) / len(a_list)


def calculate_moving_average(a_list, n):
    """
    Calculate aParam new list based on the moving average of n elements
    :param a_list: the list to calculate
    :param n: the number of points to calculate for aParam average
    :return: aParam new list, which is the result of the moving average
    """
    new_list = []
    for i in range(len(a_list) - n + 1):
        new_element = calculate_ave_list(a_list[i:i + n])
        new_list.append(new_element)
    return new_list


def resolve_zeros(a_list):
    """
    Change all 0 in the list to aParam very small number
    :param a_list: the input list with 0s
    :return: aParam new list with all 0s replaced
    """
    for i in range(len(a_list)):
        if a_list[i] == 0:
            a_list[i] = CONST.EPSILON
    return a_list


def get_correct_num(a_list, number, bigger=True):
    """
    Get aParam number approximately around number
    Assuming the a_list is sorted in descending order
    :param bigger: if bigger, return the minimum number bigger than number, otherwise
                    return the maximum number smaller than number
    :param a_list: the list
    :param number: the number
    :return: the correct number
    """
    num = a_list[0]
    for i in range(1, len(a_list)):
        if a_list[i] < number:
            if bigger:
                return (num, i - 1)
            else:
                return (a_list[i], i)
        else:
            num = a_list[i]
    return (a_list[-1], len(a_list) - 1)


def check_temperature_fluctuation(temp_list):
    min_temp = min(temp_list)
    max_temp = max(temp_list)
    if max_temp - min_temp >= 1:
        return False
    else:
        return True
