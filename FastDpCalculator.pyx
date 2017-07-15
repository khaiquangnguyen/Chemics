from libc.math cimport exp
from scipy import *
import numpy as np
cimport numpy as np
import scipy.constants
DTYPE = np.float64
ctypedef np.float64_t DTYPE_t


def cal_cc(dp, lambda_air):
    return 1 + 2 * lambda_air / dp * (1.257 + 0.4 * exp(-1.1 * dp / 2 / lambda_air))

def find_dp(long double dp, double lambda_air, int n):
    cdef long double dp_old = dp * 1 * n
    cdef long double c
    cdef long double dp_new
    cdef int i
    cdef double t
    for i in range(1000):
        t =  exp(-1.1 * dp_old / 2 / lambda_air)
        c =  1 + 2 * lambda_air / dp_old * (1.257 + 0.4 * t)
        dp_new = dp * c * n
        if abs(dp_new - dp_old) / dp_old < 0.000001:
            break
        else:
            dp_old = dp_new
    return dp_new


def calculate_fraction(np.ndarray a_list, int charge_number, np.ndarray coef_list=None):
    """
    :param a_list: the input list
    :param coef_list: the list of coefficient
    :return: the list of fraction after calculation
    """
    new_list = []
    cdef double epsilon = 0.0000001
    cdef double e = scipy.constants.e
    cdef double e0 = scipy.constants.epsilon_0
    cdef double k = scipy.constants.k
    cdef double t = scipy.constants.zero_Celsius + 25
    cdef double z = 0.875
    cdef double p = 1013
    cdef double nair = 0.000001458 * t ** 1.5 / (t + 110.4)
    cdef double lambda_air = 2 * nair / 100 / p / (8 * 28.84 / pi / 8.314 / t) ** 0.5 * 1000 ** 0.5

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
            f1 = e / sqrt(4 * pi ** 2 * e0 * i * pow(10,-9) * k * t)
            f2 = -((charge_number - (2 * pi * e0 * i * pow(10,-9) * k * t / e ** 2) * log(z)) ** 2)
            f3 = (4 * pi * e0 * i * pow(10,-9) * k * t / e ** 2)
            value = f1 * exp(f2 / f3)
            new_list.append(value)
    return np.asarray(new_list)


def correct_charges(np.ndarray[DTYPE_t] particle_diameter_list,
                    np.ndarray[DTYPE_t] cn_list,
                    np.ndarray[DTYPE_t] ccn_list,
                    np.ndarray[DTYPE_t] cn_fixed_list,
                    np.ndarray[DTYPE_t] ccn_fixed_list,
                    np.ndarray[DTYPE_t] g_cn_list,
                    np.ndarray[DTYPE_t] g_ccn_list):
# def correct_charges( particle_diameter_list,
#                      cn_list,
#                      ccn_list,
#                      cn_fixed_list,
#                     ccn_fixed_list,
#                     g_cn_list,
#                     g_ccn_list):

    cdef double asymp = 99999
    cdef double epsilon = 0.0000000001
    cdef double e = scipy.constants.e
    cdef double e0 = scipy.constants.epsilon_0
    cdef double k = scipy.constants.k
    cdef double t = scipy.constants.zero_Celsius + 25
    cdef double z = 0.875
    cdef double p = 1013
    cdef double nair = 0.000001458 * pow(t,1.5) / (t + 110.4)
    cdef double lambdaAir = 2 * nair / 100 / p / pow((8 * 28.84 / pi / 8.314 / t),0.5) * pow(1000,0.5)
    coeficientList = [[-0.0003, -0.1014, 0.3073, -0.3372, 0.1023, -0.0105],
                      [-2.3484, 0.6044, 0.48, 0.0013, -0.1553, 0.032],
                      [-44.4756, 79.3772, -62.89, 26.4492, -5.748, 0.5049]]
    cdef int i
    cdef int j
    cdef int n
    cdef double maxUpperBinBound
    cdef double moveDoubletCounts
    cdef double moveTripletCounts
    cdef double lowerBinBound
    cdef double upperBinBound
    # frac0List = calculate_fraction(diameterList,0,coeficientList[0])
    cdef np.ndarray[DTYPE_t] frac1List = calculate_fraction(particle_diameter_list, 1, np.asarray(coeficientList[1]))
    cdef np.ndarray[DTYPE_t] frac2List = calculate_fraction(particle_diameter_list, 2, np.asarray(coeficientList[2]))
    cdef np.ndarray[DTYPE_t] frac3List = calculate_fraction(particle_diameter_list, 3)
    chargeList = []
    for i in particle_diameter_list:
        aDList = [0]
        for k in range(1, 4):
            # c = cal_cc(i * pow(10,-9), lambdaAir)
            c = i * pow(10,-9)
            c = 1 + 2 * lambdaAir / c * (1.257 + 0.4 * exp(-1.1 * c / 2 / lambdaAir))

            dp = pow(10,9) * find_dp(i * pow(10,-9) / c, lambdaAir, k)
            aDList.append(dp)
        chargeList.append(aDList)
    # second part of correct charges
    cn_fixed_list = np.copy(cn_list)
    ccn_fixed_list = np.copy(ccn_list)
    maxUpperBinBound = (particle_diameter_list[-1] + particle_diameter_list[-2]) / 2
    lenDpList = len(particle_diameter_list)
    for i in range(lenDpList):
        n = lenDpList - i - 1
        moveDoubletCounts = frac2List[n] / (frac1List[n] + frac2List[n] + frac3List[n]) * cn_list[n]
        moveTripletCounts = frac3List[n] / (frac1List[n] + frac2List[n] + frac3List[n]) * cn_list[n]
        cn_fixed_list[n] = cn_fixed_list[n] - moveDoubletCounts - moveTripletCounts
        ccn_fixed_list[n] = ccn_fixed_list[n] - moveDoubletCounts - moveTripletCounts
        if chargeList[n][2] <= maxUpperBinBound:
            j = lenDpList - 2
            while (True):
                upperBinBound = (particle_diameter_list[j] + particle_diameter_list[j + 1]) / 2
                lowerBinBound = (particle_diameter_list[j] + particle_diameter_list[j - 1]) / 2
                if upperBinBound > chargeList[n][2] >= lowerBinBound:
                    cn_fixed_list[j] = cn_fixed_list[j] + moveDoubletCounts
                    if chargeList[n][2] < asymp:
                        if g_ccn_list[j] > epsilon:
                            ccn_fixed_list[j] = ccn_fixed_list[j] + moveDoubletCounts * \
                                                                              g_ccn_list[j] / \
                                                                              g_cn_list[j]
                    else:
                        ccn_fixed_list[j] = ccn_fixed_list[j] + moveDoubletCounts
                    break
                j -= 1

        if chargeList[n][3] < maxUpperBinBound:
            j = lenDpList - 2
            while (True):
                upperBinBound = (particle_diameter_list[j] + particle_diameter_list[j + 1]) / 2
                lowerBinBound = (particle_diameter_list[j] + particle_diameter_list[j - 1]) / 2
                if upperBinBound > chargeList[n][3] >= lowerBinBound:
                    cn_fixed_list[j] = cn_fixed_list[j] + moveTripletCounts
                    if chargeList[n][3] < asymp:
                        ccn_fixed_list[j] = ccn_fixed_list[j] + moveTripletCounts * ccn_list[j] / cn_list[j]
                    else:
                        ccn_fixed_list[j] = ccn_fixed_list[j] + moveTripletCounts
                    break
                j -= 1
    for i in range(len(ccn_fixed_list)):
        if ccn_fixed_list[i] / cn_fixed_list[i] < -0.01:
            ccn_fixed_list[i] = 0

    g_ccn_list = np.copy(ccn_fixed_list)
    g_cn_list = np.copy(cn_fixed_list)
    return (particle_diameter_list,cn_list,ccn_list,cn_fixed_list,ccn_fixed_list,g_cn_list,g_ccn_list)




