###############################################################################
# tvg.py
# Author: Morag Brown
###############################################################################

import IPython,random,math
from itertools import islice
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import argparse
import re
import statistics
import csv

def write_to_file(filename, write_type = 'a', write = 'row', data = ['']):
        '''Writes/appends data to specified csv file - writes an empty row if write_type, write and data aren't specified'''
        with open(filename, write_type, encoding='UTF8', newline='') as f:
            writer = csv.writer(f)
            if write == 'row':
                writer.writerow(data)
            elif write == 'rows':
                writer.writerows(data)

if __name__ == '__main__':
    # command line arg parsing
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='generates input test vector for either/both test cases')
    p.add_argument('--fs', dest='fs', type=int, action='store', help='sample rate in MSPS')
    p.add_argument('--f1', dest='f1', type=int, action='store', help='frequency of first sinusoid in MHz (max <= fs/2)')
    p.add_argument('--f2', dest='f2', type=int, action='store', help='frequency of second sinusoid in MHz (max <= fs/2)')
    p.add_argument('--tvg0', dest='tvg0', action='store_true', help='case 0 - radio astronomy')
    p.add_argument('--tvg1', dest='tvg1', action='store_true', help='case 1 - general')
    args = p.parse_args()

    plt.style.use('seaborn') # I personally prefer seaborn for the graph style, but you may choose whichever you want.
    params = {"ytick.color" : "black",
              "xtick.color" : "black",
              "axes.labelcolor" : "black",
              "axes.edgecolor" : "black",
              "text.usetex" : True,
              "font.family" : "serif",
              "font.serif" : ["Computer Modern Serif"]}
    plt.rcParams.update(params)

    input_path = "./inputs"
    calc_file = input_path + "/" + "tvg_calc_checks.csv"
    write_to_file(calc_file, 'w')

    time = 100                             # length of plots
    N = 2**16                              # length of input vector = longest FFT
    fs = args.fs*1e6                       # sample rate
    t = np.arange(0, N/fs, 1/fs)           # time vector
    freqs = [args.f1*1e6, args.f2*1e6]     # sine frequency
    
    if args.tvg0:
        write_to_file(calc_file, 'a', 'row', ["Test Vector Generator 0 Calculations"])

        noise_rms = 1/8                    # rms of noise
        mu, sigma = 0, noise_rms
        noise = np.random.normal(mu, sigma, N)
        noise_rms_calc = math.sqrt(sum(noise**2)/N)

        snr = -10                          # power snr in dB
        a_rms = noise_rms/math.sqrt(10)  # rms of sine sinusoid based on snr
        for f in freqs:
            sine = math.sqrt(2)*a_rms*np.sin(2*np.pi*f*t)
            
            sine_rms_calc = math.sqrt(sum(sine**2)/N)
            snr_actual = 20*np.log10(sine_rms_calc/noise_rms_calc)
            
            test_vector = noise + sine
            spectrum = abs(np.fft.fft(test_vector))
            
            time_energy = sum(test_vector**2)
            time_ave_power = time_energy/N
            
            freq_energy = (1/N)*sum(spectrum**2)
            freq_ave_power = (1/N)*freq_energy
            
            psd = (1/N**2)*(spectrum**2)
            normalised_psd = psd/freq_ave_power
            psd_check = sum(normalised_psd)     # should be = 1

            spect_snr = 10*np.log10(2*normalised_psd.max()/(sum(normalised_psd)-2*normalised_psd.max()))

            header = ["expected snr", "actual snr (time)", "actual snr (freq)", "noise rms", "actual noise rms", "sine rms", "actual sine rms", "time average power", "freq average power"] 
            data = [snr, round(snr_actual, 4), round(spect_snr, 4) , noise_rms, noise_rms_calc, a_rms, sine_rms_calc, time_ave_power, freq_ave_power]

            filename = "tvg0_{}msps_{}mhz_-10dBsnr".format(args.fs, round(f/(1e6)))
            write_to_file(calc_file, 'a', 'row', [filename])
            write_to_file(calc_file, 'a', 'row', header)
            write_to_file(calc_file, 'a', 'row', data)
            write_to_file(calc_file)

            # plot
            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, constrained_layout=True)
            fig.suptitle("Test Vector Generator 0")

            ax1.plot(t[0:time], noise[0:time])
            ax1.set(title = "Generated White Gaussian Noise (WGN)")
            ax2.plot(t[0:time], sine[0:time])
            ax2.set(title = "Generated {}MHz Sinusoid".format(round(f/(1e6))))
            ax2.set(ylabel = "Amplitude")
            ax3.plot(t[0:time], test_vector[0:time])
            ax3.set(title = "Test Vector = WGN + Sinusoid")
            ax3.set(xlabel = "Time")

            plt.savefig(input_path + "/" + filename + ".png")
            plt.show()
            plt.close()

            np.savetxt(input_path + "/" + filename, test_vector)

    if args.tvg1:
        write_to_file(calc_file, 'a', 'row', ["Test Vector Generator 1 Calculations"])

        lsb = 2**(-15)
        mu, sigma = 0, 0.5*lsb
        noise = np.random.normal(mu, sigma, N)
        noise_rms_calc = math.sqrt(sum(noise**2)/N)

        a = 0.9        
        for f in freqs:
            sine = 0.9*np.sin(2*np.pi*f*t)

            sine_rms_calc = math.sqrt(sum(sine**2)/N)
            snr_calc = 20*np.log10(sine_rms_calc/noise_rms_calc)

            test_vector = sine + noise
            spectrum = abs(np.fft.fft(test_vector))
            
            time_energy = sum(test_vector**2)
            time_ave_power = time_energy/N
            
            freq_energy = (1/N)*sum(spectrum**2)
            freq_ave_power = (1/N)*freq_energy
            
            psd = (1/N**2)*(spectrum**2)
            normalised_psd = psd/freq_ave_power
            psd_check = sum(normalised_psd)     # should be = 1

            spect_snr = 10*np.log10(2*normalised_psd.max()/(sum(normalised_psd)-2*normalised_psd.max()))

            header = ["snr (time)", "snr (freq)", "noise rms", "actual noise rms", "sine max amplitude", "sine rms", "time average power", "freq average power"]
            data = [round(snr_calc, 4), round(spect_snr, 4) , sigma, noise_rms_calc, a, round(sine_rms_calc, 4), time_ave_power, freq_ave_power]

            filename = "tvg1_{}msps_{}mhz_{}dBsnr".format(args.fs, round(f/(1e6)), round(snr_calc))
            write_to_file(calc_file, 'a', 'row', [filename])
            write_to_file(calc_file, 'a', 'row', header)
            write_to_file(calc_file, 'a', 'row', data)
            write_to_file(calc_file)

            # plot
            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, constrained_layout=True)
            fig.suptitle("Test Vector Generator 1")

            ax1.plot(t[0:time], noise[0:time])
            ax1.set(title = "Generated White Gaussian Noise (WGN)")
            ax2.plot(t[0:time], sine[0:time])
            ax2.set(title = "Generated {}MHz Sinusoid".format(round(f/(1e6))))
            ax2.set(ylabel = "Amplitude")
            ax3.plot(t[0:time], test_vector[0:time])
            ax3.set(title = "Test Vector = WGN + Sinusoid")
            ax3.set(xlabel = "Time")

            plt.savefig(input_path + "/" + filename + ".png")
            plt.show()
            plt.close()
 
            np.savetxt(input_path + "/" + filename, test_vector)


















