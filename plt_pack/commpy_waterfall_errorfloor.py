__author__ = 'yihanjiang'


import matplotlib.pylab as plt

if __name__ == '__main__':
    SNRS = [-1.    ,     -0.71346414 ,-0.41715244 ,-0.11037421 , 0.20763712,  0.53773555,
            0.88087633,  1.23813262 , 1.61071576,  2.   ,  2.45850188,  2.94256209   ]
    BER_AWGN = [2.11690000e-03  , 3.31300000e-04  , 7.24000000e-05 ,  2.80000000e-05,
           1.31000000e-05  , 6.50000000e-06 ,  3.40000000e-06  , 1.60000000e-06,
           8.00000000e-07 ,  6.00000000e-07 ,1.80000000e-07  , 8.00000000e-08]

    plt.figure(1)
    plt.yscale('log')
    plt.plot(SNRS, BER_AWGN)
    plt.show()