__author__ = 'yihanjiang'
import matplotlib.pylab as plt
if __name__ == '__main__':
    Turbo_bursty =  [0.0292, 0.0226, 0.0158, 0.0302, 0.033, 0.0452, 0.0418, 0.0356, 0.04, 0.0744, 0.0668, 0.0954, 0.068, 0.0332, 0.0354, 0.041, 0.0318, 0.0434, 0.0428, 0.0426, 0.0414, 0.077, 0.0404, 0.0694, 0.0404, 0.0382, 0.0136, 0.0388, 0.0924, 0.0466, 0.0412, 0.067, 0.0448, 0.0394, 0.0926, 0.1092, 0.1576, 0.053, 0.064, 0.0682, 0.0694, 0.059, 0.0626, 0.0536, 0.1114, 0.0768, 0.1054, 0.1578, 0.1116, 0.128, 0.4922, 0.1352, 0.2008, 0.0846, 0.077, 0.044, 0.0582, 0.071, 0.0848, 0.0742, 0.041, 0.0356, 0.0376, 0.043, 0.065, 0.056, 0.0488, 0.0696, 0.047, 0.04, 0.0826, 0.0318, 0.0352, 0.0356, 0.0288, 0.0278, 0.0354, 0.0322, 0.0264, 0.0392, 0.037, 0.0356, 0.0388, 0.0836, 0.0468, 0.0426, 0.0158, 0.0916, 0.075, 0.0446, 0.0334, 0.0478, 0.0392, 0.0298, 0.0372, 0.0408, 0.0718, 0.0572, 0.0602, 0.0468]
    Turbo_non_bursty  = [0.0032, 0.005, 0.0028, 0.0058, 0.0074, 0.0064, 0.008, 0.009, 0.008, 0.007, 0.0082, 0.007, 0.0078, 0.0082, 0.0072, 0.008, 0.0054, 0.0046, 0.0078, 0.0082, 0.0076, 0.0062, 0.0072, 0.006, 0.007, 0.0062, 0.0018, 0.005, 0.0082, 0.0066, 0.008, 0.0084, 0.0078, 0.007, 0.0084, 0.009, 0.0072, 0.0058, 0.007, 0.0082, 0.013, 0.013, 0.0104, 0.0064, 0.0152, 0.0076, 0.0074, 0.0132, 0.0064, 0.0058, 0.0054, 0.005, 0.006, 0.006, 0.0064, 0.0044, 0.0082, 0.0052, 0.0058, 0.0066, 0.0056, 0.007, 0.0068, 0.0078, 0.0114, 0.0082, 0.0062, 0.0118, 0.0084, 0.0074, 0.0068, 0.0064, 0.0084, 0.006, 0.006, 0.0046, 0.0082, 0.0068, 0.0068, 0.0088, 0.0074, 0.005, 0.0068, 0.008, 0.0102, 0.01, 0.0032, 0.0078, 0.0086, 0.0096, 0.007, 0.0066, 0.0086, 0.0074, 0.0098, 0.006, 0.008, 0.0096, 0.013, 0.0096]

    plt.figure(1)
    plt.title('Turbo Decoding Error Propagation on Block Length 100')
    plt.xlabel('Code Block Position')
    plt.ylabel('Positional BER')
    p1, = plt.plot(Turbo_bursty, 'r',label = 'Positional BER Radar Bursty Noise')
    p2, = plt.plot(Turbo_non_bursty, label = 'Positional BER Radar AWGN')
    plt.legend(handles = [p1, p2])
    plt.show()