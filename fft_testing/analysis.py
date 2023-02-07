###############################################################################
# analysis.py
# Author: Morag Brown
###############################################################################

import IPython
import os
from os import listdir
from os.path import isfile, join
import csv 
from itertools import islice
import matplotlib.pyplot as plt 
import argparse
import re
import statistics
from scipy.stats import pearsonr
import numpy as np
from subprocess import call

def write_to_file(results_file, write_type = 'a', write = 'row', data = ['']):
        '''Writes/appends data to specified csv file - writes an empty row if write_type, write and data aren't specified'''
        with open(results_file, write_type, encoding='UTF8', newline='') as f:
            writer = csv.writer(f)
            if write == 'row':
                writer.writerow(data)
            elif write == 'rows':
                writer.writerows(data)

def load_fftw_output(fftw_outputs_path, input_file, nfft):
        '''Loads FFTW output for current input file'''
        fftw_output_file = join(fftw_outputs_path, 'fftw_' + str(nfft) + "_" + input_file + '.output')

        with open(fftw_output_file) as f:
            lines_in = f.readlines()
        fftw_output = np.array([line.split('\n')[0] for line in lines_in], dtype=float)

        return fftw_output


def load_fft_output(test_type, outputs_path, hw_type, fft_type, nfft, input_file):
        all_output_files = [f for f in listdir(outputs_path) if isfile(join(outputs_path, f))]

        output_files = [f for f in all_output_files if (hw_type + "_" + fft_type + "_" + str(nfft) + "_" + input_file) in f]
        outputs = list()

        for output_file in output_files:
            m = re.search('sr-(.+?).output', output_file)
            sr = m.group(1)
            with open(join(outputs_path, output_file)) as f:
                lines_in = f.readlines()

            if test_type == 'sim':
                tmp = list([line.split('\n')[0] for line in lines_in])
                output = list([output_file.replace('.output', ''), 'shift: {}'.format(sr), tmp[0], tmp[1], tmp[2], np.array(tmp[3:], dtype=float)])
                outputs.append(output)

            elif test_type == 'hw':
                tmp = list([line.split('\n')[0] for line in lines_in])
                output = list([output_file.replace('.output', ''), 'shift: {}'.format(sr), np.array(tmp, dtype=float)])
                outputs.append(output)

        return outputs

def calc_power(nfft, scale, spectrum, psd_file, plot_title):
        '''Calculates average power, normalised PSD and SPDR of spectrum + saves and plots normalised PSD'''
        # IPython.embed()
        pxx = (1/nfft)*sum(spectrum**2)
        ave_pxx = (1/nfft)*pxx
        ave_scaled_pxx = ave_pxx*scale**2

        psd = (1/nfft)*(spectrum**2)
        normalised_psd = psd/pxx

        tmp = np.sort(normalised_psd)[::-1]
        sfdr = 10*np.log10(tmp[0]) - 10*np.log10(tmp[1])

        # snr for only pos frequencies components
        spect_snr = 10*np.log10(normalised_psd.max()/(sum(normalised_psd)-normalised_psd.max()))

        freq_axis = np.arange(0, nfft//2)*(800/nfft)
 
        fig, (ax1, ax2) = plt.subplots(2, sharex=True)
        fig.suptitle("Normalised PSD: {}".format(plot_title))
        ax1.plot(freq_axis, normalised_psd, 'cornflowerblue')
        ax1.set(ylabel = "$P_{xx}[f] \ (linear)$")
        ax2.plot(freq_axis, 10*np.log10(normalised_psd), 'cornflowerblue')
        ax2.set(ylabel = "$P_{xx}[f] \ (dB)$")
        ax2.set(xlabel = "Frequency (MHz)")

        plt.savefig(psd_file + ".png")
        np.savetxt(psd_file, normalised_psd)
        plt.close()

        return ave_pxx, ave_scaled_pxx, psd, normalised_psd, sfdr, spect_snr


def calc_mse(nfft, fft_psd, fftw_psd):
        '''Calculates mean squared error of FFT and FFTW (normalised, ideally) PSDs'''
        squared_diff = (fft_psd - fftw_psd)**2
        mse = (1/nfft)*sum(squared_diff)

        return mse


def sim_analysis(fftw_output, outputs_path, input_file, nfft, fftw_results_file, results_file, shift_dict, demux):
        test_type = 'sim'

        # fftw stuff
        fftw_header = ["Average Power", "Average Power (dB)", "SFDR (dB)", "SNR (dB)"] 
        fftw_psd_file = join(os.path.split(fftw_results_file)[0], "normalised_psd_fftw_nfft-" + str(nfft))
        fftw_plot_title = "FFTW" + "\n" + "FFT Length = {}".format(nfft)
        scale = 1
        fftw_ave_power, _, fftw_psd, fftw_normalised_psd, fftw_sfdr, fftw_snr = calc_power(nfft, scale, fftw_output, fftw_psd_file, fftw_plot_title)
        write_to_file(fftw_results_file, 'a', 'row', fftw_header)
        write_to_file(fftw_results_file, 'a', 'row', [fftw_ave_power, 10*np.log10(fftw_ave_power), fftw_sfdr, fftw_snr])
        write_to_file(fftw_results_file)

        hw_types = ['au50', 'skarab']
        au50_data = list()
        skarab_data = list()
        for hw_type in hw_types:
           if hw_type == 'au50':
               fft_types = ['astron', 'casper', 'xilinx']
               for fft_type in fft_types:
                   outputs = load_fft_output(test_type, outputs_path, hw_type, fft_type, nfft, input_file)
                   results = list()
                   for output in outputs:
                       fft_psd_file = join(os.path.split(results_file)[0], "{}in_normalised_psds/".format(demux) + output[0] + "_{}in_normalised_psd".format(demux))
                       m = re.search('sr-(.+?)_', fft_psd_file)
                       sr = m.group(1)
                       if fft_type == 'astron':
                           scale = 4*shift_dict[nfft][sr]
                       else:
                           scale = shift_dict[nfft][sr]
                       fft_plot_title = "[" + hw_type.upper() + "] " + fft_type.upper() + "\n" + "FFT Length = {}".format(nfft) + ", Shifting Schedule = {}".format(sr) + ", L = {}".format(demux)
                       fft_ave_power, fft_scaled_ave_power, fft_psd, fft_normalised_psd, fft_sfdr, fft_snr = calc_power(nfft, scale, output[5], fft_psd_file, fft_plot_title)
                       fft_mse = calc_mse(nfft, fft_normalised_psd, fftw_normalised_psd)
                       sim_time, latency, overflow = output[2].split(": ")[1], output[3].split(": ")[1], output[4].split(": ")[1]
                       data = [fft_type, sr, fft_ave_power, fft_scaled_ave_power, 10*np.log10(fft_scaled_ave_power), fft_mse, fft_sfdr, fft_snr, float(sim_time), (latency), (overflow)]
                       au50_data.append(data)

           elif hw_type == 'skarab':
               fft_types = ['astron', 'casper']
               for fft_type in fft_types:
                   outputs = load_fft_output(test_type, outputs_path, hw_type, fft_type, nfft, input_file)
                   results = list()
                   for output in outputs:
                       fft_psd_file = join(os.path.split(results_file)[0], "{}in_normalised_psds/".format(demux) + output[0] + "_{}in_normalised_psd".format(demux))
                       m = re.search('sr-(.+?)_', fft_psd_file)
                       sr = m.group(1)
                       if fft_type == 'astron':
                           scale = 4*shift_dict[nfft][sr]
                       else:
                           scale = shift_dict[nfft][sr]
                       fft_plot_title = "[" + hw_type.upper() + "] " + fft_type.upper() + "\n" + "FFT Length = {}".format(nfft) + ", Shifting Schedule = {}".format(sr) + ", L = {}".format(demux)
                       fft_ave_power, fft_scaled_ave_power, fft_psd, fft_normalised_psd, fft_sfdr, fft_snr = calc_power(nfft, scale, output[5], fft_psd_file, fft_plot_title)
                       fft_mse = calc_mse(nfft, fft_normalised_psd, fftw_normalised_psd)
                       sim_time, latency, overflow = output[2].split(": ")[1], output[3].split(": ")[1], output[4].split(": ")[1]
                       data = [fft_type, sr, fft_ave_power, fft_scaled_ave_power, 10*np.log10(fft_scaled_ave_power), fft_mse, fft_sfdr, fft_snr, float(sim_time), (latency), (overflow)]
                       skarab_data.append(data)

        # write results to file
        header = ["FFT Core", "Shifting Schedule", "Average Power", "Scaled Average Power", "SAP dB", "MSE", "SFDR (dB)", "SNR", "Sim Time", "Latency", "Overflow?"]
        write_to_file(results_file, 'a', 'row', ["AU50"])
        write_to_file(results_file, 'a', 'row', header)
        write_to_file(results_file, 'a', 'rows', au50_data)
        write_to_file(results_file, 'a', 'row', ["SKARAB"])
        write_to_file(results_file, 'a', 'row', header)
        write_to_file(results_file, 'a', 'rows', skarab_data) 
        write_to_file(results_file)


def hw_analysis(fftw_output, outputs_path, input_file, nfft, fftw_results_file, results_file, shift_dict, demux):
        test_type = 'hw'

        fftw_psd_file = join(os.path.split(fftw_results_file)[0], "normalised_psd_fftw_nfft-" + str(nfft))
        fftw_plot_title = "FFTW" + "\n" + "FFT Length = {}".format(nfft)
        scale = 1
        fftw_ave_power, _, fftw_psd, fftw_normalised_psd, fftw_sfdr, fftw_snr = calc_power(nfft, scale, fftw_output, fftw_psd_file, fftw_plot_title)

        hw_types = ['au50', 'skarab']
        au50_data = list()
        skarab_data = list()
        for hw_type in hw_types:
           if hw_type == 'au50':
               fft_types = ['astron', 'casper', 'xilinx']
               for fft_type in fft_types:
                   outputs = load_fft_output(test_type, outputs_path, hw_type, fft_type, nfft, input_file)
                   results = list()
                   for output in outputs:
                       fft_psd_file = join(os.path.split(results_file)[0], "{}in_normalised_psds/".format(demux) + output[0] + "_{}in_normalised_psd".format(demux))
                       m = re.search('sr-(.+?)_', fft_psd_file)
                       sr = m.group(1)
                       if fft_type == 'astron':
                           scale = 4*shift_dict[nfft][sr]
                       else:
                           scale = shift_dict[nfft][sr]
                       fft_plot_title = "[" + hw_type.upper() + "] " + fft_type.upper() + "\n" + "FFT Length = {}".format(nfft) + ", Shifting Schedule = {}".format(sr) + ", L = {}".format(demux)
                       fft_ave_power, fft_scaled_ave_power, fft_psd, fft_normalised_psd, fft_sfdr, fft_snr = calc_power(nfft, scale, output[2], fft_psd_file, fft_plot_title)
                       fft_mse = calc_mse(nfft, fft_normalised_psd, fftw_normalised_psd)
                       data = [fft_type, sr, fft_ave_power, fft_scaled_ave_power, 10*np.log10(fft_scaled_ave_power), fft_mse, fft_sfdr, fft_snr]
                       au50_data.append(data)

           elif hw_type == 'skarab':
               fft_types = ['astron', 'casper']
               for fft_type in fft_types:
                   outputs = load_fft_output(test_type, outputs_path, hw_type, fft_type, nfft, input_file)
                   results = list()
                   for output in outputs:
                       fft_psd_file = join(os.path.split(results_file)[0], "{}in_normalised_psds/".format(demux) + output[0] + "_{}in_normalised_psd".format(demux))
                       m = re.search('sr-(.+?)_', fft_psd_file)
                       sr = m.group(1)
                       if fft_type == 'astron':
                           scale = 4*shift_dict[nfft][sr]
                       else:
                           scale = shift_dict[nfft][sr]
                       fft_plot_title = "[" + hw_type.upper() + "] " + fft_type.upper() + "\n" + "FFT Length = {}".format(nfft) + ", Shifting Schedule = {}".format(sr) + ", L = {}".format(demux)
                       fft_ave_power, fft_scaled_ave_power, fft_psd, fft_normalised_psd, fft_sfdr, fft_snr = calc_power(nfft, scale, output[2], fft_psd_file, fft_plot_title)
                       fft_mse = calc_mse(nfft, fft_normalised_psd, fftw_normalised_psd)
                       data = [fft_type, sr, fft_ave_power, fft_scaled_ave_power, 10*np.log10(fft_scaled_ave_power), fft_mse, fft_sfdr, fft_snr]
                       skarab_data.append(data)

        # write results to file
        header = ["FFT Core", "Shifting Schedule", "Average Power", "Scaled Average Power", "SAP dB", "MSE", "SFDR (dB)", "SNR"]
        write_to_file(results_file, 'a', 'row', ["AU50"])
        write_to_file(results_file, 'a', 'row', header)
        write_to_file(results_file, 'a', 'rows', au50_data)
        write_to_file(results_file, 'a', 'row', ["SKARAB"])
        write_to_file(results_file, 'a', 'row', header)
        write_to_file(results_file, 'a', 'rows', skarab_data)
        write_to_file(results_file)

# START
if __name__ == '__main__':
    # command line arg parsing
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='does the thing')
    p.add_argument('--inputs', dest='inputs', type=str, action='store', help='path to input files')
    p.add_argument('--sim', dest='sim', action='store_true', help='analyses simulation outputs if specified')
    p.add_argument('--hw', dest='hw', action='store_true', help='analyses hardware outputs if specified')
    args = p.parse_args()

    plt.style.use('seaborn') 
    params = {"ytick.color" : "black",
              "xtick.color" : "black",
              "axes.labelcolor" : "black",
              "axes.edgecolor" : "black",
              "text.usetex" : True,
              "font.family" : "serif",
              "font.serif" : ["Computer Modern Serif"]}
    plt.rcParams.update(params)

    N = 100                                        # length of time vector for input plot
    fftw_outputs_path = './fftw/outputs'

    inputs_path = args.inputs
    input_files = [f for f in listdir(inputs_path) if isfile(join(inputs_path, f))]

    for input_file in input_files:
        m = re.search('_(.+?)msps', input_file)
        if m:
            fs = int(m.group(1))*1e6 

        # create results directory for test vector if it doesn't already exist (vector generator should do this probably)
        results_path = "./results/" + input_file
        if not os.path.exists(results_path):
            os.makedirs(results_path)

        # create separate directory for fftw results (normalised PSDs + csv file containing ave power, SFDR for each nfft)
        fftw_results_path = results_path + "/fftw"
        if not os.path.exists(fftw_results_path):
            os.makedirs(fftw_results_path)
        fftw_results_file = join(fftw_results_path, "fftw_results.csv")
        write_to_file(fftw_results_file, 'w', 'row', ["Input Signal: {}".format(input_file)])
        write_to_file(fftw_results_file)
        
        # iterate through outputs for both demux factors
        demux_vals = [4, 8]
        # demux_vals = [4]
        for demux in demux_vals:

            # create directories and results files for specified test type
            if args.sim:
                sim_outputs_path = "./simulation/outputs/{}in".format(demux)
                sim_results_path = results_path + "/{}in/simulation".format(demux)
                sim_plots_path = join(sim_results_path, "{}in_normalised_psds".format(demux))
                if not os.path.exists(sim_results_path):
                    os.makedirs(sim_results_path)
                if not os.path.exists(sim_plots_path):
                    os.makedirs(sim_plots_path)
                sim_results_file = join(sim_results_path, "{}in_simulation_results.csv".format(demux))
                write_to_file(sim_results_file, 'w', 'row', ["Input Signal: {}".format(input_file)])
                write_to_file(sim_results_file)

            if args.hw:
                hw_outputs_path = "./hardware/outputs/{}in".format(demux)
                hw_results_path = results_path + "/{}in/hardware".format(demux) 
                hw_plots_path = join(hw_results_path, "{}in_normalised_psds".format(demux))
                if not os.path.exists(hw_results_path):
                    os.makedirs(hw_results_path)
                if not os.path.exists(hw_plots_path):
                    os.makedirs(hw_plots_path)
                hw_results_file = join(hw_results_path, "{}in_hardware_results.csv".format(demux))
                write_to_file(hw_results_file, 'w', 'row', ["Input Signal: {}".format(input_file)])
                write_to_file(hw_results_file)

            # scale val is how many times shifted [eg 1k alt is 1010101010 = 682 in dec, but only 5 downshifts so scale by 32]
            # no shifting gets scale factor of 1
            shift_dict = {2**10: {'none':1, 'alt':32, 'all':1023}, 2**14: {'none':1, 'alt':128, 'all':16383}, 2**16: {'none':1, 'alt':256, 'all':65535}}

            fft_lengths = [2**10, 2**14, 2**16]
            # fft_lengths = [2**10]
            for nfft in fft_lengths:
                write_to_file(fftw_results_file, 'a', 'row', ["FFT Length = {}".format(nfft)])
                # load fftw reference
                fftw_output = load_fftw_output(fftw_outputs_path, input_file, nfft)

                if args.sim:
                     write_to_file(sim_results_file, 'a', 'row', ["FFT Length = {}".format(nfft)])
                     sim_analysis(fftw_output, sim_outputs_path, input_file, nfft, fftw_results_file, sim_results_file, shift_dict, demux)
 
                if args.hw:
                     write_to_file(hw_results_file, 'a', 'row', ["FFT Length = {}".format(nfft)])
                     hw_analysis(fftw_output, hw_outputs_path, input_file, nfft, fftw_results_file, hw_results_file, shift_dict, demux)
