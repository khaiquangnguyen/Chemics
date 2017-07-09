
cdef extern from "vfastexp.h":
    double exp_approx "EXP" (double)


cdef double cal_cc(double dp, double lambda_air)except? -2:
    cdef double t =  exp_approx(-1.1 * dp / 2 / lambda_air)
    cdef double m =  1 + 2 * lambda_air / dp * (1.257 + 0.4 * t)
    return m

def find_dp(double dp, double lambda_air, int n):
    cdef double dp_old = dp * 1 * n
    cdef double c
    cdef double dp_new
    cdef int i
    for i in range(1000):
        c = cal_cc(dp_old, lambda_air)
        dp_new = dp * c * n
        if abs(dp_new - dp_old) / dp_old < 0.000001:
            break
        else:
            dp_old = dp_new
    return dp_new
