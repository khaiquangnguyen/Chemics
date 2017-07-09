
cdef extern from "vfastexp.h":
    double exp_approx "EXP" (double)

def find_dp(double dp, double lambda_air, int n):
    cdef double dp_old = dp * 1 * n
    cdef double c
    cdef double dp_new
    cdef int i
    cdef double t
    for i in range(1000):
        t =  exp_approx(-1.1 * dp_old / 2 / lambda_air)
        c =  1 + 2 * lambda_air / dp * (1.257 + 0.4 * t)
        dp_new = dp * c * n
        if abs(dp_new - dp_old) / dp_old < 0.000001:
            break
        else:
            dp_old = dp_new
    return dp_new
