import scipy.optimize as opt
import scipy.stats
import scipy.signal
from HelperFunctions import *
import FastDpCalculator
import settings as CONST

class Scan:
    """
    Class for one scan, which contains all the information related to a scan
    """

    def __init__(self):
        # status of the scan. 1 for good, 0 for not
        self.status = 1
        # status code. The reason why the scan is not good
        self.status_code = 0
        # the flow rate
        self.counts_to_conc = 0.0
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
        # diameter midpoints.
        self.diameter_midpoints = []
        # ave diameter from smps file
        self.ave_smps_diameters = []
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
        self.true_super_sat = 0
        self.processed_ave_ccnc_sizes = []
        self.processed_normalized_concs = []
        # smps and ccnc data after we correct for the charges
        self.corrected_smps_counts = []
        self.corrected_ccnc_counts = []
        # necessary parameters for sigmoid fit
        # begin rise. end rise, begin asymp, end asymp
        # the main reason we need this is because we have to handle fitting in two sigmoid lines
        self.sigmoid_params = []
        # b,d and c
        self.functions_params = []
        # dps of interest
        # dp50, dp50wet, dp50+20
        self.dps = []
        # sigmoid line y values
        self.sigmoid_y_vals = []

    def pre_align_self_test(self):
        # if the length of smps data is not in sync with duration, the consider the scan invalid
        if len(self.raw_smps_counts) != self.duration:
            self.status = 0
            self.set_status_code(2)
            # todo: perform Hartigan's dip test to test for bimodality

    def post_align_self_test(self):
        # todo: implement something so that we know that a scan is good or not
        # check for error in super saturation
        for i in range(len(self.processed_super_sats)):
            if not are_floats_equal(self.true_super_sat, self.processed_super_sats[i]):
                self.true_super_sat = None
                self.set_status(0)
                self.set_status_code(7)
                break
        if numpy.std(self.processed_T1s) > 1.5 or  numpy.std(self.processed_T2s) > 1.5 or  numpy.std(
                self.processed_T3s) > 1.5:
                self.set_status(0)
                self.set_status_code(7)

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
        self.raw_smps_counts = numpy.asarray(self.raw_smps_counts, dtype=float64)
        self.raw_ccnc_counts = numpy.asarray(self.raw_ccnc_counts, dtype=float64)
        self.raw_normalized_concs = numpy.asarray(self.raw_normalized_concs)
        self.diameter_midpoints = numpy.asarray(self.diameter_midpoints)
        self.raw_ave_ccnc_sizes = numpy.asarray(self.raw_ave_ccnc_sizes)
        self.raw_super_sats = numpy.asarray(self.raw_super_sats)
        self.ave_smps_diameters = numpy.asarray(self.ave_smps_diameters)

    def align_smps_ccnc_data(self, base_shift_factor):
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
        # got to normalized the dndlogdp list
        self.processed_normalized_concs = normalize_dndlogdp_list(self.raw_normalized_concs)
        # transform them first
        ccnc_counts = self.raw_ccnc_counts
        t1s = self.raw_T1s
        t2s = self.raw_T2s
        t3s = self.raw_T3s
        super_sats = self.raw_super_sats
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
            super_sats = super_sats[self.shift_factor:]
            ave_ccnc_sizes = ave_ccnc_sizes[self.shift_factor:]
        # if shift factor is negative
        else:
            # populate ccnc counts with 0s in the fronts
            ccnc_counts = fill_zeros_to_begin(ccnc_counts, abs(self.shift_factor))
            t1s = fill_zeros_to_begin(t1s, abs(self.shift_factor))
            t2s = fill_zeros_to_begin(t2s, abs(self.shift_factor))
            t3s = fill_zeros_to_begin(t3s, abs(self.shift_factor))
            super_sats = fill_zeros_to_begin(super_sats, abs(self.shift_factor))
            ave_ccnc_sizes = fill_zeros_to_begin(ave_ccnc_sizes, abs(self.shift_factor))
        # If we still don't have enough data, populate with 0s
        ccnc_counts = fill_zeros_to_end(ccnc_counts, self.duration)
        t1s = fill_zeros_to_end(t1s, self.duration)
        t2s = fill_zeros_to_end(t2s, self.duration)
        t3s = fill_zeros_to_end(t3s, self.duration)
        super_sats = fill_zeros_to_end(super_sats, self.duration)
        ave_ccnc_sizes = fill_zeros_to_end(ave_ccnc_sizes, self.duration)
        self.processed_ccnc_counts = ccnc_counts[:self.duration]
        self.processed_T1s = t1s[:self.duration]
        self.processed_T2s = t2s[:self.duration]
        self.processed_T3s = t3s[:self.duration]
        self.processed_super_sats = super_sats[:self.duration]
        self.true_super_sat = self.processed_super_sats[0]
        self.processed_ave_ccnc_sizes = ave_ccnc_sizes[:self.duration]
        # perform self test
        self.post_align_self_test()

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
            ave_smps_dp, smps, ccnc, corrected_smps, corrected_ccnc, prev_smps, prev_ccnc = \
                FastDpCalculator.correct_charges(
                ave_smps_dp, smps, ccnc, corrected_smps, corrected_ccnc, prev_smps, prev_ccnc)
        # save the smps and ccnc data after charge correction to new variables
        self.corrected_smps_counts = corrected_smps
        self.corrected_ccnc_counts = corrected_ccnc

    def helper_get_ratio_corrected_smooth(self,asym_limits = [0.75,1.5]):
        ratio_corrected = []
        ccnc = self.corrected_ccnc_counts
        smps = self.corrected_smps_counts
        # smooth the corrections
        ccnc = heavy_smooth(ccnc)
        smps = heavy_smooth(smps)
        in_the_beginning = True
        # we take care of the case where there can be unreasonable ccnc data at the beginning
        for i in range(len(self.corrected_smps_counts)):
            if smps[i] > max(smps) // 20:
                in_the_beginning = False
            if in_the_beginning:
                ratio_corrected.append(0)
            else:
                ratio_corrected.append(safe_div(ccnc[i], smps[i]))
        ratio_corrected = ratio_corrected[:len(self.ave_smps_diameters)]
        # smooth the hell out of the data
        ratio_corrected = heavy_smooth(ratio_corrected)
        # remove huge data points
        for i in range(len(ratio_corrected)):
            if ratio_corrected[i] > asym_limits[1]:
                if i > 0:
                    ratio_corrected[i] = ratio_corrected[i - 1]
                else:
                    ratio_corrected[i] = 0
        return ratio_corrected

    def cal_params_for_sigmoid_fit(self):
        asym_limits = [0.75, 1.5]
        while True:
            try:
                self.cal_params_for_sigmoid_fit_single_loop(asym_limits)
                break
            except:
                # widen the gap between the limits so that we can cover a larger area
                asym_limits = [asym_limits[0]-0.1, asym_limits[1] + 0.1]

    def cal_params_for_sigmoid_fit_single_loop(self, asym_limits = [0.75,1.5]):
        if not self.is_valid():
            return
        # first, got to clean up everything
        self.sigmoid_params = []
        self.functions_params = []
        self.dps = []
        self.sigmoid_y_vals = []
        # this is just an approximation. It doesn't have to be really accurate
        ratio_corrected = self.helper_get_ratio_corrected_smooth(asym_limits)
        # next, find the min_dp by finding the first point at which all points before it are less than 0.1,
        # with an error of 5%
        # set up all the params we need
        begin_rise = 0
        # get the bottom inflation point. Great!
        for i in range(len(ratio_corrected) - 1, 1, -1):
            value = ratio_corrected[i]
            pre_values = ratio_corrected[:i]
            pre_percent = cal_percentage_less(pre_values, value, 0.05)
            if value < 0.15 and pre_percent > 0.80:
                begin_rise = self.ave_smps_diameters[i]
                break
        # get the top inflation point
        # get the top 10% of the data
        high_index_list = []
        for i in range(len(ratio_corrected)):
            if asym_limits[0] <= ratio_corrected[i] <= asym_limits[1]:
                high_index_list.append(i)
        high_list = ratio_corrected[high_index_list]
        high_index_list = outliers_iqr_ver_2(high_list, high_index_list)
        # if the return is an error/meaning we can't detect any outliers
        # then we simply get the highest possible value from the array
        if high_index_list == -1:
            high_index_list = numpy.argmax(ratio_corrected)
        high_index_list = numpy.sort(high_index_list)
        # next, remove outliers
        top_inflation_index = high_index_list[0]
        last_dp_index = high_index_list[-1]
        # now, chose the last dp so that the std of the remaining data is within 0.1
        while True:
            if len(high_index_list) <= 2:
                break
            last_dp_index = high_index_list[-1]
            interested_list = ratio_corrected[top_inflation_index:last_dp_index]
            std = numpy.std(interested_list)
            if std <= (asym_limits[1] - asym_limits[0])/2:
                break
            else:
                high_index_list = high_index_list[:-1]
        end_rise = self.ave_smps_diameters[top_inflation_index]
        end_asymp = self.ave_smps_diameters[last_dp_index]
        self.sigmoid_params.append([begin_rise,end_rise,end_rise, end_asymp])

    def fit_sigmoids(self, asym_limits = None):
        if not self.is_valid():
            return
        self.functions_params = []
        self.sigmoid_y_vals = []
        self.dps = []
        for i in range(len(self.sigmoid_params)):
            self.fit_one_sigmoid(i)

    def fit_one_sigmoid(self, params_set_index, asym_limits = [0.75, 1.5]):
        begin_rise = self.sigmoid_params[params_set_index][0]
        end_rise = self.sigmoid_params[params_set_index][1]
        begin_asymp = self.sigmoid_params[params_set_index][2]
        end_asymp = self.sigmoid_params[params_set_index][3]
        ratio_corrected = safe_div_array(self.corrected_ccnc_counts, self.corrected_smps_counts)
        # smooth the hell out of the data
        ratio_corrected = heavy_smooth(ratio_corrected)
        # get b
        ave_list = []
        for i in range(len(self.ave_smps_diameters)):
            if begin_asymp < self.ave_smps_diameters[i] < end_asymp and asym_limits[0] < ratio_corrected[i] < asym_limits[1]:
                ave_list.append(ratio_corrected[i])
        b = get_ave_none_zero(ave_list)
        if not asym_limits[0] <= b <= asym_limits[1]:
            b = 1
        def fn(x, d, c):
            return b / (1 + (x / d) ** c)
        x_list = []
        y_list = []
        # get all data points on the rise
        for i in range(len(self.ave_smps_diameters)):
            if begin_rise < self.ave_smps_diameters[i] < end_rise:
                x_list.append(self.ave_smps_diameters[i])
                y_list.append(ratio_corrected[i])
        # get all data points on the asymp
        for i in range(len(self.ave_smps_diameters)):
            if begin_asymp < self.ave_smps_diameters[i] < end_asymp:
                x_list.append(self.ave_smps_diameters[i])
                y_list.append(ratio_corrected[i])
        x_list = numpy.asarray(x_list)
        y_list = numpy.asarray(y_list)
        try:
            result = opt.curve_fit(fn, x_list, y_list, bounds=([begin_rise, -200], [end_asymp+1, -1]), method="trf")
            d = result[0][0]
            c = result[0][1]
        except:
            d = 60
            c = -2
        self.functions_params.append([b,d,c])
        dp_50 = d
        dp_50_wet = 0
        dp_50_20_wet = 0
        # find dp50 and dp50 wet
        for i in range(1,len(self.ave_smps_diameters)):
            if self.ave_smps_diameters[i] > d:
                dp_50_wet = self.processed_ave_ccnc_sizes[i-1]
                break
        # find dp50 +20 (wet)
        for i in range(1, len(self.ave_smps_diameters)):
            if self.ave_smps_diameters[i] > d + 20:
                dp_50_20_wet = self.processed_ave_ccnc_sizes[i - 1]
                break
        sigmoid_points = [0]
        for i in range(1,len(self.ave_smps_diameters)):
            if begin_rise <= self.ave_smps_diameters[i] <= end_asymp:
                sigmoid_points.append(fn(self.ave_smps_diameters[i],d,c))
            else:
                sigmoid_points.append(sigmoid_points[i-1])
        self.sigmoid_y_vals.append(sigmoid_points)
        dp_50 = round(dp_50,3)
        dp_50_wet = round(dp_50_wet,3)
        dp_50_20_wet = round(dp_50_20_wet,3)
        self.dps.append([dp_50,dp_50_wet,dp_50_20_wet])

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

    def decode_status_code(self):
        if self.status_code == 0:
            return "no error"
        elif self.status_code == 1:
            return "There is no equivalent ccnc data for this scan. This is likely because ccnc data start at a later" \
                   "time than smps scan, or it ends before the smps scan."
        elif self.status_code == 2:
            return "The length of SMPS for this scan does not agree with the indicated scan duration of the experiment."
        elif self.status_code == 3:
            return "The distribution of SMPS and CCNC data has a low correlation with those of the next few scans."
        elif self.status_code == 4:
            return "The program can not locate the reference point for SMPS data."
        elif self.status_code == 5:
            return "The program can not locate the reference point for CCNC data."
        elif self.status_code == 6:
            return "The scan do not have enough CCNC or SMPS data. Most likely because we shift the data too much"
        elif self.status_code == 7:
            return "The super saturation rate do not remain constant throughout the scan!"
        elif self.status_code == 8:
            return "The temperature do not remain constant enough throughout the scan!"
        elif self.status_code == 9:
            return "The Scan is manually disabled by the user!"

            # return "The super saturation rate does not remain constant throughout the scan duration."

    def add_to_raw_normalized_concs(self, new_data):
        self.raw_normalized_concs.append(float(new_data))

    def add_to_raw_smps_counts(self, new_data):
        self.raw_smps_counts.append(int(new_data) * self.counts_to_conc)

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

    def add_to_diameter_midpoints(self, new_data):
        self.diameter_midpoints.append(float(new_data))

    def set_index_in_ccnc_data(self, index):
        self.index_in_ccnc_data = index

    def add_to_ave_smps_diameters(self, new_data):
        self.ave_smps_diameters.append(new_data)

    def set_status(self, status):
        self.status = status

    def set_status_code(self, code):
        self.status_code = code

    def set_shift_factor(self, factor):
        if factor < len(self.raw_ccnc_counts):
            self.shift_factor = factor
        else:
            self.shift_factor = 0

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

    def set_counts_2_conc(self, n):
        self.counts_to_conc = n

    def set_sigmoid_params(self,params):
        self.sigmoid_params = params
        self.fit_sigmoids()
