import numpy as np
import StringIO

#
# CASCADES.DAT
#
__cascades = """
North cascades data from Hosking.

Site ID        N    Mean    L-CV  L-skew   L-kur      T5
--------------------------------------------------------

350304        98  19.685  0.1209  0.0488  0.1433 -0.0004
351433        59  62.580  0.0915  0.0105  0.1569  0.0020
351862        90  40.852  0.1124  0.0614  0.1541 -0.0058
351897        61  46.045  0.1032  0.0417  0.1429 -0.0022
352997        65  45.021  0.0967 -0.0134  0.1568  0.0173
353445        86  31.042  0.1328 -0.0176  0.1206  0.0235
353770        78  80.143  0.1008  0.0943  0.1967  0.0856
356907        72  41.305  0.1143  0.0555  0.1210  0.0487
357169        67  30.585  0.1107  0.0478  0.1371  0.0316
357331        99  32.932  0.1179  0.0492  0.0900  0.0225
357354        49  17.560  0.1308  0.0940  0.1273  0.0352
358466        61  69.518  0.1119 -0.0429  0.0927 -0.0061
450945        69  47.653  0.1018  0.0435  0.1446 -0.0056
451233        73 102.501  0.1025  0.0182  0.1047 -0.0221
453284        70  52.413  0.1054 -0.0224  0.1664  0.0035
454764        66  79.696  0.1174  0.0124  0.1317 -0.0176
454769        59  44.643  0.1115 -0.0346  0.1032  0.0083
457773        74  58.655  0.1003  0.0446  0.1450 -0.0379
458773        82  39.024  0.1046  0.0128  0.1583  0.0443
"""

dt = [('id', 'a6'), ('N', 'i'), ('lmom', 'f', (5,))]
input = np.loadtxt(StringIO.StringIO(__cascades), dtype=dt, skiprows=5)

#
# CASCADES.OUT
#

__results_1 = """
    North Cascades                                         19 SITES

 SITE    N      NAME       L-CV   L-SKEW  L-KURT   D(I)
    1   98  350304        0.1209  0.0488  0.1433   0.60
    2   59  351433        0.0915  0.0105  0.1569   1.02
    3   90  351862        0.1124  0.0614  0.1541   0.38
    4   61  351897        0.1032  0.0417  0.1429   0.23
    5   65  352997        0.0967 -0.0134  0.1568   0.93
    6   86  353445        0.1328 -0.0176  0.1206   2.63
    7   78  353770        0.1008  0.0943  0.1967   2.12
    8   72  356907        0.1143  0.0555  0.1210   0.45
    9   67  357169        0.1107  0.0478  0.1371   0.11
   10   99  357331        0.1179  0.0492  0.0900   1.61
   11   49  357354        0.1308  0.0940  0.1273   2.08
   12   61  358466        0.1119 -0.0429  0.0927   1.52
   13   69  450945        0.1018  0.0435  0.1446   0.31
   14   73  451233        0.1025  0.0182  0.1047   1.30
   15   70  453284        0.1054 -0.0224  0.1664   1.58
   16   66  454764        0.1174  0.0124  0.1317   0.29
   17   59  454769        0.1115 -0.0346  0.1032   1.04
   18   74  457773        0.1003  0.0446  0.1450   0.43
   19   82  458773        0.1046  0.0128  0.1583   0.38
"""

discordance = np.loadtxt(StringIO.StringIO(__results_1), skiprows=4, usecols=[6])

"""

     WEIGHTED MEANS       0.1103  0.0279  0.1366

 PARAMETERS OF REGIONAL KAPPA DISTRIBUTION   0.9542  0.1533  0.1236 -0.2955


 ***** HETEROGENEITY MEASURES *****
 (NUMBER OF SIMULATIONS  =   500)

 OBSERVED     S.D. OF GROUP L-CV          =  0.0104
 SIM. MEAN OF S.D. OF GROUP L-CV          =  0.0095
 SIM. S.D. OF S.D. OF GROUP L-CV          =  0.0016
 STANDARDIZED TEST VALUE H(1)             =  0.62

 OBSERVED AVE.  OF L-CV / L-SKEW DISTANCE =  0.0339
 SIM. MEAN OF AVE. L-CV / L-SKEW DISTANCE =  0.0447
 SIM. S.D. OF AVE. L-CV / L-SKEW DISTANCE =  0.0072
 STANDARDIZED TEST VALUE H(2)             = -1.49

 OBSERVED AVE.  OF L-SKEW/L-KURT DISTANCE =  0.0405
 SIM. MEAN OF AVE. L-SKEW/L-KURT DISTANCE =  0.0578
 SIM. S.D. OF AVE. L-SKEW/L-KURT DISTANCE =  0.0073
 STANDARDIZED TEST VALUE H(3)             = -2.37


 ***** GOODNESS-OF-FIT MEASURES *****
 (NUMBER OF SIMULATIONS  =   500)

 GEN. LOGISTIC        L-KURTOSIS= 0.167   Z VALUE=  3.46
 GEN. EXTREME VALUE   L-KURTOSIS= 0.111   Z VALUE= -2.94
 GEN. NORMAL          L-KURTOSIS= 0.123   Z VALUE= -1.51 *
 PEARSON TYPE III     L-KURTOSIS= 0.123   Z VALUE= -1.60 *
 GEN. PARETO          L-KURTOSIS= 0.006   Z VALUE=-14.75


 PARAMETER ESTIMATES FOR DISTRIBUTIONS ACCEPTED AT THE 90% LEVEL

 GEN. NORMAL          0.994  0.195 -0.057
 PEARSON TYPE III     1.000  0.196  0.171
 WAKEBY               0.568  2.003  7.330  0.244 -0.270

 QUANTILE ESTIMATES
                      0.010  0.020  0.050  0.100  0.200  0.500  0.900  0.950  0.990  0.999
 GEN. NORMAL          0.569  0.616  0.688  0.753  0.834  0.994  1.254  1.331  1.480  1.654
 PEARSON TYPE III     0.570  0.616  0.688  0.753  0.834  0.994  1.254  1.331  1.480  1.653
 WAKEBY               0.590  0.610  0.666  0.740  0.840  0.993  1.259  1.341  1.483  1.603
"""
rmom = [0.1103, 0.0279, 0.1366]
vobs = [0.0104, .0339, .0405]
vbar = [.0095, .0447, .0578]
vsd = [.0016, .0072, .0073]
h = [.62, -1.49, -2.37]
z = [3.46, -2.94, -1.51, -1.60, -14.75]
para_gen_normal = [.994, .195, -.057]
para_pearson = [1., .196, .171]
para_wakeby = [0.568, 2.003, 7.330, 0.244, -0.270]
