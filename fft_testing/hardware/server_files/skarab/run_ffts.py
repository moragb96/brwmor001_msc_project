import IPython,casperfpga,time,struct,os,sys,logging,pylab,random,math,re,argparse
import numpy as np
import matplotlib.pyplot as plt
from itertools import islice
from os import listdir
from os.path import isfile, join

def get_input_files(inputs_path, fpgs_path):
        # returns lists of input and fpg files based on paths passed as command line arguments
        logging.info("Getting list of fpg and input files") 

        fpg_files = [f for f in listdir(fpgs_path) if isfile(join(fpgs_path, f))]
        input_files = [f for f in listdir(inputs_path) if isfile(join(inputs_path, f))]

        demux_dict = {'4in': 4, '8in': 8}
        demux = demux_dict[fpgs_path.split('/')[2]]

        return fpg_files, input_files, demux


def parse_fpg_file(fpg_file):
        # determines FFT implementation type and FFT length based on fpg file name
        logging.info("Determining FFT length and type from fpg file")

        fft_type = fpg_file.split('_')[0]
        fft_length = fpg_file.split('_')[4]

        length_dict = {'1k' : 2**10, '16k' : 2**14, '64k' : 2**16}
        nfft = length_dict[fft_length]

        return fft_type, nfft


def parse_input_file(inputs_path, input_file, fft_type, nfft, sr, sr_type, demux):
        # reads signal from input file and generates output file namei
        logging.info("Loading test vector from file")

        with open(join(inputs_path,input_file)) as f1: 
            lines_in = f1.readlines()
            input_sig = np.array([line.split()[0] for line in lines_in])
    
        sig = input_sig.astype(float)[0:nfft]

        m = re.search('_(.+?)msps', input_file)
        if m:
            fs = int(m.group(1))*1e6

        scaled_sig = np.empty(nfft)
        for i in range(len(sig)):
           scaled_sig[i] = sig[i]*(2**15)

        output_file = "./outputs/{}in/skarab_{}_{}_".format(str(demux), fft_type, nfft) + input_file + "_sr-{}".format(sr_type) + ".output"

        return sig, fs, scaled_sig.astype(int), output_file


def write_signal(sig, demux, n):
        # create dictionary of demuxed arrays 
        logging.info("Writing input signal to input BRAMs")

        d = {}
        for i in range(demux):
             d['d{}'.format(i)] = sig[i:n:demux]    

        # create dictionary of packed buffers to be written to BRAM and write
        buf = {}
        for j in range(demux):
            buf['buf{}'.format(j)] = struct.pack('>{}i'.format(str(n/demux)), *d['d{}'.format(j)])    # n/demux -> number of samples to each BRAM
            fpga.write('bram_in{}'.format(j), buf['buf{}'.format(j)], 0)

        return d, buf


def read_input_brams(demux, n, plot = True):
        logging.info("Reading from input BRAMs")

        b = {}
        for i in range(demux):
            b['in{}'.format(i)] = np.array(struct.unpack('>{}i'.format(str(n/demux)), fpga.read('bram_in{}'.format(i) , (n/demux)*4 , 0)), dtype = np.float)

        interleave_in = []
        for j in range(n/demux):
            for jj in range(demux):
                interleave_in.append(b['in{}'.format(jj)][j])

        if plot:
            logging.info("Plotting input BRAMs")
            plt.plot(interleave_in)
            plt.show()

        return np.array(interleave_in, dtype=np.float)

def read_input_ss(demux, plot = True):
        # snapshots capturing data from input BRAMs are same dimensions for all FFT lengths (256 deep, 16*demux wide)
        logging.info("Reading from input snapshot")
        
        ss_in = fpga.snapshots.snapshot_input_ss.read(arm=False)['data']
        ss_sig = []
        for i in range(256):
            for ii in range(demux):
                ss_sig.append(ss_in['in{}'.format(ii)][i])

        if plot:
            logging.info("Plotting input snapshot")
            plt.plot(ss_sig)
            plt.show()

        return np.array(ss_sig, dtype=np.float)

def get_output_data(demux, n, fft, output_file, fs, plot = True):
        # reads output of FFT from snapshots
        logging.info("Reading output snapshots")
        
        re_out = fpga.snapshots.snapshot_re_out_ss.read(arm=False)['data']
        im_out = fpga.snapshots.snapshot_im_out_ss.read(arm=False)['data']

        interleave_re = []
        interleave_im = []
        # casper
        if fft == 'casper':
            for i in range(n/demux):
                for ii in range(demux/2):
                    interleave_re.append(re_out['re{}'.format(ii)][i])
                    interleave_im.append(im_out['im{}'.format(ii)][i])

            interleave_re_array = np.array(interleave_re,dtype=np.float)
            interleave_im_array = np.array(interleave_im,dtype=np.float)

        # astron
        # astron outputs only pos freq channels, but as re0 = [a0, b0, a1, b1 ...] etc
        if fft == 'astron':
            if demux == 4:
                for j in range(demux):
                    # with delay into snapshots, need to index from 1 to get b vals
                    interleave_re.append(re_out['re{}'.format(j)][1:(n/demux):2])
                    interleave_im.append(im_out['im{}'.format(j)][1:(n/demux):2])

            elif demux == 8:
                re_out1 = fpga.snapshots.snapshot_re_out1_ss.read(arm=False)['data']
                im_out1 = fpga.snapshots.snapshot_im_out1_ss.read(arm=False)['data']

                for j in range(demux/2):
                    interleave_re.append(re_out['re{}'.format(j)][1:(n/demux):2])
                    interleave_im.append(im_out['im{}'.format(j)][1:(n/demux):2])
                for j in range(demux/2, demux):
                    interleave_re.append(re_out1['re{}'.format(j)][1:(n/demux):2])
                    interleave_im.append(im_out1['im{}'.format(j)][1:(n/demux):2])

            interleave_re_array = np.array(interleave_re,dtype=np.float)
            interleave_re_array = interleave_re_array.flatten()
            interleave_im_array = np.array(interleave_im,dtype=np.float)
            interleave_im_array = interleave_im_array.flatten() 

        complex_vals = interleave_re_array + interleave_im_array*1j
        complex_out = abs(interleave_re_array + interleave_im_array*1j)

        if plot:
            f_axis = np.linspace(0,nfft/2 - 1,nfft/2)*fs/nfft
            plt.plot(f_axis, complex_out[0:nfft/2])
            plt.title("SKARAB: {} FFT, {} points".format(fft.upper(), str(nfft)))
            plt.xlabel("Freq (MHz)")
            plt.show()

        np.savetxt(output_file + ".complex", complex_vals)
        np.savetxt(output_file, complex_out)
        #IPython.embed()

        return complex_out

# START
if __name__ == '__main__':
        
        # command line arguments
        p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='runs fft tests on hardware and captures results to file')
        p.add_argument('--board', dest='board', type=str, action='store', default='', help='hardware platform hostname or port')
        p.add_argument('--inputs', dest='inputs', type=str, action='store', help='path to input files')
        p.add_argument('--fpgs', dest='fpgs', type=str, action='store', help='path to fpg files')
        args = p.parse_args()

        logging.basicConfig(level=eval('logging.INFO'))

        # connect to board
        logging.info('Connecting to SKARAB')
        fpga = casperfpga.CasperFpga(args.board)
        time.sleep(0.2)

        if fpga.is_connected():
            logging.info('Connected to SKARAB')
        else:
            logging.info('ERROR connecting to SKARAB')
            exit_fail()

        # load list of fpg and input files, and extract demux value
        fpg_files, input_files, demux = get_input_files(args.inputs, args.fpgs)

        # run each fft for every input file (tvg signal is length of longest fft)
        for input_file in input_files:
            logging.info('---------------------------------------------')
            logging.info('Using input file: {}'.format(input_file))

            # iterate over fpg files for each fft (and fft length)
            for fpg_file in fpg_files:
                logging.info('---------------------------------------------')
                logging.info('Using fpg file: {}'.format(fpg_file))
                fft_type, nfft = parse_fpg_file(fpg_file)
    
                # program board or load fpg information 
                logging.info('---------------------------------------------')
                logging.info('Programming SKARAB')
                fpga.upload_to_ram_and_program(join(args.fpgs, fpg_file))
                logging.info('Done programming SKARAB')
    
                #  shift_type = ['none', 'lastH', 'firstH', 'alt', 'all']
                # shift_dict = {2**10: [0, 31, 992, 682, 1023], 2**14: [0, 127, 16256, 10922, 16383], 2**16: [0, 255, 65280, 43690, 65535]}
                shift_type = ['none', 'alt']
                shift_dict = {2**10: [0, 682], 2**14: [0, 10922], 2**16: [0, 43690]}
                #shift_type = ['all']
                #shift_dict = {2**10: [1023], 2**14: [16383], 2**16: [65535]}
                shift_schedules = shift_dict[nfft]

                for isr, sr in enumerate(shift_schedules):
                    sr_type = shift_type[isr]
                    logging.info('---------------------------------------------')
                    logging.info('Using shifting schedule: {}'.format(sr_type))
    
                    # parse input file according to nfft (only take nfft worth of tvg signal to write to bram)
                    sig, fs, scaled_sig, output_file = parse_input_file(args.inputs, input_file, fft_type, nfft, sr, sr_type, demux)
                    logging.info('Output will be written to: {}'.format(output_file))    
                    logging.info('---------------------------------------------')

                    #IPython.embed()
    
                    # configure design for each input file 
                    logging.info('Configuring design')
    
                    logging.info('Arming snapshots')
                    fpga.snapshots.snapshot_input_ss.arm()
                    fpga.snapshots.snapshot_re_out_ss.arm()
                    fpga.snapshots.snapshot_im_out_ss.arm()
                    if demux == 8 and fft_type == 'astron':
                         fpga.snapshots.snapshot_re_out1_ss.arm()
                         fpga.snapshots.snapshot_im_out1_ss.arm()
    
                    d, buf = write_signal(scaled_sig, demux, nfft)

                    fpga.write_int('in_enable', 0)
                    fpga.write_int('reset', 1)
                    fpga.write_int('reset', 0)
    
                    fpga.write_int('shift', sr)
                    fpga.write_int('in_enable', 1)
    
                    #read_input_brams(demux=4, n=nfft, plot = True)
                    #read_input_ss(demux=4, plot = True)
    
                    fft_out = get_output_data(demux, nfft, fft_type, output_file, fs, plot = True)
    
                    fpga.write_int('in_enable', 0)
                    fpga.write_int('reset', 1)
                    fpga.write_int('reset', 0)
                    time.sleep(2)    

                    logging.info('---------------------------------------------')
    
                logging.info('---------------------------------------------')
#            IPython.embed()
