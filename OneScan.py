from MainView import *
import matplotlib.pyplot as plt
import matplotlib
import scipy.optimize as opt
import scipy.stats
import scipy.signal
from Exceptions import *
import gc
from HelperFunctions import *
import FastDpCalculator
import settings as CONST
import matplotlib.patches as mpatches
import os
from PySide.QtCore import SIGNAL
from PySide.QtCore import *

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
matplotlib.style.use('ggplot')


class OneScan:
    """
    Class for one scan, which contains all the information related to a scan
    """

    def __init__(self):
        # status of the scan. 1 for good, 0 for not
        self.status = 1
        # status code. The reason why the scan is not good
        self.status_code = 0
        # the flow rate
        self.counts_to_conc = 0.06
        # the scan #
        self.index = 0
        # the position of the scan in ccnc data
        self.index_in_ccnc_data = 0
        # what time the scan starts. Format is hh:mm:ss
        self.start_time = None
        # what time the scan ends. Format is hh:mm:ss
        self.end_time = None
        # duration of the scan
        self.duration = 0
        # up and down scan time. Very useful to align the data
        self.scan_up_time = 0
        self.scan_down_time = 0
        # the list of super saturation rate.
        self.raw_super_sats = []
        # the temperature measured. Called "T1 Read, T2 Read and T3 Read" in the CCNC file
        self.raw_T1s = []
        self.raw_T2s = []
        self.raw_T3s = []
        # SMPS and CCNC counts. SMPS count is calculated from smps file and ccnc count is calculated from ccnc file
        self.raw_smps_counts = []
        self.raw_ccnc_counts = []
        # calculated average size of ccnc particles from ccnc file
        self.raw_ave_ccnc_sizes = []
        # normalized concentration. Also called dN/dlogDp.
        self.raw_normalized_concs = []
        # The reference point of smps and ccnc data. Normally, this is the local minimum between the peaks
        # Used to align smps and ccnc data
        self.ref_index_smps = 0
        self.ref_index_ccnc = 0
        self.shift_factor = 0
        # Processed data. Since we only need a portion of the raw data and work with it, we don't use all the raw data
        # Furthermore, we need to perform certain processing techniques on the data
        self.processed_smps_counts = []
        self.processed_ccnc_counts = []
        self.processed_T1s = []
        self.processed_T2s = []
        self.processed_T3s = []
        self.processed_super_sats = []
        self.processed_ave_ccnc_sizes = []
        self.processed_normalized_concs = []
        # smps and ccnc data after we correct for the charges
        self.corrected_smps_counts = []
        self.corrected_ccnc_counts = []
        self.ave_smps_diameters = []

    def pre_align_self_test(self):
        # if the length of smps data is not in sync with duration, the consider the scan invalid
        if len(self.raw_smps_counts) != self.duration:
            self.status = 0
            self.set_status_code(2)
            # todo: perform Hartigan's dip test to test for bimodality

    def post_align_self_test(self):
        # todo: implement something so that we know that a scan is good or not
        # check for temperature problem
        # x_axis = numpy.arange(len(self.raw_T1s))
        # scipy.stats.linregress(x_axis, self.raw_T1s)
        # x_axis = numpy.arange(len(self.raw_T2s))
        # scipy.stats.linregress(x_axis, self.raw_T2s)
        # x_axis = numpy.arange(len(self.raw_T3s))
        # scipy.stats.linregress(x_axis, self.raw_T3s)
        # check for super saturation problem
        return 1

    def do_basic_trans(self):
        """
        Perform basic transformation before we try to align the smps and ccnc data.
        We want to keep the raw data intact so that we can use them later when we need them
        :return:
        """
        # convert everything to numpy arrays
        # got to generate processed_smps_counts early

        self.raw_T1s = numpy.asarray(self.raw_T1s)
        self.raw_T2s = numpy.asarray(self.raw_T2s)
        self.raw_T3s = numpy.asarray(self.raw_T3s)
        self.raw_smps_counts = numpy.asarray(self.raw_smps_counts,dtype=float64)
        self.raw_ccnc_counts = numpy.asarray(self.raw_ccnc_counts, dtype=float64)
        self.raw_normalized_concs = numpy.asarray(self.raw_normalized_concs)
        self.raw_ave_ccnc_sizes = numpy.asarray(self.raw_ave_ccnc_sizes)
        self.raw_super_sats = numpy.asarray(self.raw_super_sats)
        self.ave_smps_diameters = numpy.asarray(self.ave_smps_diameters)

    def align_smps_ccnc_data(self,base_shift_factor):
        # the assumption here is that the shift factor is always positive
        # the base shift factor is used a reference. The general idea is that the shift factors of two scans should
    # be within 1 index of each other
        self.shift_factor = base_shift_factor
        ccnc_counts = self.raw_ccnc_counts
        smps_counts = self.raw_smps_counts
        # if we have an invalid scan, then do nothing.
        if not self.is_valid():
            return base_shift_factor
        # normalize ccnc counts. Method is based on what is stated
        if len(self.raw_ccnc_counts > 3):
            ccnc_counts = smooth(ccnc_counts)
        # try to find the ref index in smps data. This is an easy one, almost always accurate
        self.ref_index_smps = find_ref_index_smps(smps_counts, self.scan_up_time)
        if self.ref_index_smps is None:
            self.status = 0
            self.set_status_code(4)
            return base_shift_factor
        # try to find the ref index in ccnc data. Very difficult
        self.ref_index_ccnc = find_ref_index_ccnc(ccnc_counts, self.ref_index_smps + base_shift_factor)
        # if we did not manage to find and point that fits our target
        if self.ref_index_ccnc is None:
            self.status = 0
            self.set_status_code(5)
            return base_shift_factor
        else:
            self.shift_factor = self.ref_index_ccnc - self.ref_index_smps
        return self.shift_factor

    def generate_processed_data(self):
        # process smps is straightforward
        # process ccnc datas is difficult, since it also based on shift factor
        # The general idea is that, if we have negative shift factor, insert 0s before the ccnc data
        # Else, we shift the ccnc data by shift_factor points
        # Then, if we don't have enough data, we will just populate the data with 0s
        # Then we collect the data
        self.processed_smps_counts = self.raw_smps_counts
        # transform them first
        ccnc_counts = self.raw_ccnc_counts
        t1s = self.raw_T1s
        t2s = self.raw_T2s
        t3s = self.raw_T3s
        ave_ccnc_sizes = self.raw_ave_ccnc_sizes
        # if shift factor is positive
        if self.shift_factor >= 0:
            # if not enough ccnc counts to even shift, the scan is invalid
            if len(ccnc_counts) < self.shift_factor:
                self.set_status(0)
                self.set_status_code(6)
            # Get data
            ccnc_counts = ccnc_counts[self.shift_factor:]
            t1s = t1s[self.shift_factor:]
            t2s = t2s[self.shift_factor:]
            t3s = t3s[self.shift_factor:]
            ave_ccnc_sizes = ave_ccnc_sizes[self.shift_factor:]
        # if shift factor is negative
        else:
            # populate ccnc counts with 0s in the fronts
            ccnc_counts = fill_zeros_to_begin(ccnc_counts, abs(self.shift_factor))
            t1s = numpy.concatenate((numpy.asarray([0] * abs(self.shift_factor)),t1s))
            t2s = numpy.concatenate((numpy.asarray([0] * abs(self.shift_factor)), t2s))
            t3s = numpy.concatenate((numpy.asarray([0] * abs(self.shift_factor)), t3s))
            ave_ccnc_sizes =numpy.concatenate((numpy.asarray([0] * abs(self.shift_factor)), ave_ccnc_sizes))
        #  If we still don't have enough data, populate with 0s
        ccnc_counts = fill_zeros_to_end(ccnc_counts,self.duration)
        t1s = fill_zeros_to_end(t1s,self.duration)
        t2s = fill_zeros_to_end(t2s,self.duration)
        t3s = fill_zeros_to_end(t2s,self.duration)
        ave_ccnc_sizes = fill_zeros_to_end(ave_ccnc_sizes,self.duration)
        self.processed_ccnc_counts = ccnc_counts[:self.duration]
        self.processed_T1s = t1s[:self.duration]
        self.processed_T2s = t2s[:self.duration]
        self.processed_T3s = t3s[:self.duration]
        self.processed_ave_ccnc_sizes = ave_ccnc_sizes[:self.duration]

    def correct_charges(self):
        # scan is not usable, do nothing
        if not self.is_valid():
            return -1
        # otherwise, let's process it
        # first, we need to initiate some necessary variables
        ccnc = self.processed_ccnc_counts
        smps = self.processed_smps_counts
        ave_smps_dp = self.ave_smps_diameters
        # do some processing on those
        ccnc = resolve_zeros(ccnc)
        smps = resolve_zeros(smps)
        ccnc = resolve_small_ccnc_vals(ccnc)
        prev_ccnc = numpy.copy(ccnc)
        prev_smps = numpy.copy(smps)
        corrected_ccnc = numpy.copy(ccnc)
        corrected_smps = numpy.copy(smps)
        # now, correct the charges
        for i in range(CONST.NUM_OF_CHARGES_CORR):
            ave_smps_dp,smps,ccnc,corrected_smps,corrected_ccnc,prev_smps,prev_ccnc = FastDpCalculator.correct_charges(ave_smps_dp,smps,ccnc,corrected_smps,corrected_ccnc,prev_smps,prev_ccnc)
        # save the smps and ccnc data after charge correction to new variables
        self.corrected_smps_counts = smps
        self.corrected_ccnc_counts = ccnc

    def detect_bisigmoid(self):
        if not self.is_valid():
            return -1
        # basically, detect if the data may fit using two sigmoid lines
        # we know that the bisigmoid thing is caused by ccnc data having weird peaks before the dip
        # so we will do it by detecting the weirdness of ccnc data
        # extract the data first. We don't want to mess with the original data
        ccnc_data = self.processed_ccnc_counts
        # next, we know that the dip is near the scan up time index.
        # so we cut out the rest of the data
        ccnc_data = ccnc_data[:self.scan_up_time]
        # next, we cut out the very small ccnc data
        i = 0
        max_ccnc = numpy.amax(ccnc_data)
        while ccnc_data[i] < max_ccnc /10 and i < (len(ccnc_data)-1):
            i+= 1
        ccnc_data = ccnc_data[i:]
        # now, we got a decent dataset to analyze for two bumps
        # First, let's do aggressive filtering
        ccnc = self.corrected_ccnc_counts
        smps = self.processed_smps_counts
        # calculate the ratio
        ratio = []
        for i in range(len(ccnc)):
            ratio.append(safe_div(ccnc[i], smps[i]))
        ratio = numpy.asarray(scipy.signal.savgol_filter(ratio, 15, 2))
        # peaks = peakutils.indexes (ratio, min_dist = 4)
        # plt.plot(self.ave_smps_diameters,ratio)
        print ()


    def get_params_for_sigmoid_fit(self):
        return 1

    def fit_sigmoid_line(self):
        return 1

    def is_valid(self):
        return self.status == 1

    def compare_smps(self, another_scan):
        # compare two scans using correlation coefficient
        corr_coef_p = scipy.stats.pearsonr(self.raw_smps_counts, another_scan.raw_smps_counts)[1]
        # compare two scans using Kolmogorov-Smirnov statistic
        ks_stat_p = scipy.stats.ks_2samp(self.raw_smps_counts, another_scan.raw_smps_counts)[1]
        if corr_coef_p < 0.05 < ks_stat_p:
            return True
        else:
            return False

    def set_status(self, status):
        self.status = status

    def set_status_code(self, code):
        self.status_code = code

    def set_shift_factor(self,factor):
        if factor < len(self.raw_ccnc_counts):
            self.shift_factor = factor
        else:
            self.shift_factor = 0

    def decode_status_code(self):
        if self.status_code == 0:
            return "no error"
        elif self.status_code == 1:
            return "There is no equivalent ccnc data for this scan. This is likely because ccnc data start at a later" \
                   "time than smps scan, or it ends before the smps scan."
        elif self.status_code == 2:
            return "The length of SMPS for this scan does not agree with the indicated scan duration of the experiment."
        elif self.status_code == 3:
            return "The distribution of SMPS and CCNC data for this scan has a low correlation with those of the next " \
                   "" \
                   "few scans."
        elif self.status_code == 4:
            return "The program can not locate the reference point for SMPS data."
        elif self.status_code == 5:
            return "The program can not locate the reference point for CCNC data."
        elif self.status_code == 6:
            return "The scan do not have enough CCNC or SMPS data. Most likely because we shift the data too much"

            # return "The super saturation rate does not remain constant throughout the scan duration."

    def set_start_time(self, start_time):
        self.start_time = start_time

    def set_end_time(self, end_time):
        self.end_time = end_time

    def set_up_time(self, up_time):
        self.scan_up_time = up_time

    def set_down_time(self, down_time):
        self.scan_down_time = down_time

    def set_duration(self, duration):
        self.duration = duration

    def add_to_raw_normalized_concs(self, new_data):
        self.raw_normalized_concs.append(float(new_data))

    def add_to_raw_smps_counts(self, new_data):
        self.raw_smps_counts.append(int(new_data)*self.counts_to_conc)

    def add_to_raw_t1s_t2s_t3s(self, t1, t2, t3):
        t1 = float(t1)
        t2 = float(t2)
        t3 = float(t3)
        self.raw_T1s.append(t1)
        self.raw_T2s.append(t2)
        self.raw_T3s.append(t3)

    def add_to_raw_super_sats(self, ss):
        ss = float(ss)
        self.raw_super_sats.append(ss)

    def add_to_raw_ccnc_counts(self, new_data):
        self.raw_ccnc_counts.append(float(new_data))

    def add_to_raw_ave_ccnc_sizes(self, new_data):
        self.raw_ave_ccnc_sizes.append(float(new_data))

    def set_index_in_ccnc_data(self, index):
        self.index_in_ccnc_data = index

    def add_to_ave_smps_diameters(self, new_data):
        self.ave_smps_diameters.append(new_data)