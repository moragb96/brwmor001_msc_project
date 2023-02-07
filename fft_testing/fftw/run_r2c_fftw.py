###############################################################################
# run_r2c_fftw.py
# Author: Morag Brown
###############################################################################

import IPython,logging
from itertools import islice
import matplotlib.pyplot as plt
import argparse
import re
import statistics
from scipy.stats import pearsonr
import numpy as np
import time
from subprocess import call
from os import listdir
from os.path import isfile, join
import os

def parse_input_file(input_file, nfft):
        logging.info("Generating output filename")
        # generate FFTW output file names
        fftw_output_file = "./outputs/fftw_"+ str(nfft) + "_" + input_file + ".output"    # properly format output file to be passed to analysis script
        complex_output_file = fftw_output_file + ".complex"

        m = re.search('_(.+?)msps', input_file)
        if m:
            fs = int(m.group(1))*1e6

        return fftw_output_file, complex_output_file, fs

def plot_output(fftw_output_file, nfft, fs, plot):
        logging.info("Plotting FFTW output")
        with open(fftw_output_file) as f2:
            lines_out = f2.readlines()

        fftw_c_out = np.array([line.split('\n')[0] for line in lines_out], dtype=float)

        if plot:
            freq_axis = np.linspace(0, nfft//2, nfft//2)*fs/nfft # create frequency axis about 0Hz

            plt.plot(freq_axis/1e6, fftw_c_out)
            plt.title('FFTW output')
            plt.xlabel('frequency (MHz)')
        
            plt.show()

        return fftw_c_out


if __name__ == '__main__':
    # command line arg parsing
    # passed to C program that calculates FFT using FFTW
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='runs the r2c fftw executable for all given inputs and all fft lengths')
    p.add_argument('--inputs', dest='inputs', type=str, action='store', help='path to input files')
    args = p.parse_args()

    logging.basicConfig(level=eval('logging.INFO'))

    input_files = [f for f in listdir(args.inputs) if isfile(join(args.inputs, f))]

    for input_file in input_files:
        logging.info("---------------------------------------------")
        logging.info("Using input file: {}".format(input_file))

        fft_lengths = [2**10, 2**14, 2**16]
        for nfft in fft_lengths:    
            fftw_output_file, complex_output_file, fs = parse_input_file(input_file, nfft)
         
        #------------------------------FFTW COMP------------------------------#
            logging.info("Running FFTW")
            # call C program that calculates FFT of signal contained in input file using FFTW
            call(["./calc_r2c_fftw.exe", join(args.inputs, input_file), fftw_output_file, complex_output_file, str(nfft)])

            # plotting for debugging purposes
#           fftw_out = plot_output(fftw_output_file, nfft, fs, plot=True)
    
