%% AU50 (work/repos/mlib_devel): wideband factor = 4

fs = 800e6;             % sampling frequency
NFFT = 2^10;            % FFT Size
binbw = fs/NFFT *1e-6;  % BW in MHz in a single bin
N = NFFT*2;            % length of time vector
t = (1:N)/fs;           % time vector
a_sig = 0.45;         % tone amplitude
f_sig = 50*1e6;         % tone frequency 

real_sig = a_sig*cos(2*pi*f_sig*t);

sim_len = N;
SSR = 4;
bram_depth = N/SSR;
bram_width = log2(bram_depth);

% demux shared_bram
d0 = int32(real_sig(1:4:end)*(2^15));
d1 = int32(real_sig(2:4:end)*(2^15));
d2 = int32(real_sig(3:4:end)*(2^15));
d3 = int32(real_sig(4:4:end)*(2^15));

freqs = [600, 550, 500, 450, 400];

for i = 1:5
    f = string(freqs(i));
    
    filename = sprintf('/media/morag/linux_storage/storage_home/university/test_system/fft_testing/hardware/logs/au50_4in_%smhz_%s.log', f, datestr(now,'yyyy_mm_dd_HH_MM'))
    fileID = fopen(filename, 'w')
    fprintf(fileID, sprintf('%s: Starting %sMHz builds for AU50 (4 inputs)\n\n', datestr(now,'yyyy_mm_dd_HH_MM'), f))
    
    tStart = tic;

    %******************** CASPER ********************
    load_system('casper_hw_test_au50_1k')
    set_param('casper_hw_test_au50_1k/AU50', 'clk_rate', f)
    save_system('casper_hw_test_au50_1k', sprintf('./fft_testing/hardware/casper/4in/casper_hw_test_au50_1k_%smhz', f))
    set_param(sprintf('casper_hw_test_au50_1k_%smhz', f),'SimulationCommand','Update')
    jasper(sprintf('casper_hw_test_au50_1k_%smhz', f))
    save_system(sprintf('casper_hw_test_au50_1k_%smhz', f))
    close_system(sprintf('casper_hw_test_au50_1k_%smhz', f))
    
    disp(sprintf('casper_hw_test_au50_1k_%smhz DONE', f))
    disp('Waiting for 5 minutes')
    pause(300)

    load_system('casper_hw_test_au50_16k')
    set_param('casper_hw_test_au50_16k/AU50', 'clk_rate', f)
    save_system('casper_hw_test_au50_16k', sprintf('./fft_testing/hardware/casper/4in/casper_hw_test_au50_16k_%smhz', f))
    set_param(sprintf('casper_hw_test_au50_16k_%smhz', f),'SimulationCommand','Update')
    jasper(sprintf('casper_hw_test_au50_16k_%smhz', f))
    save_system(sprintf('casper_hw_test_au50_16k_%smhz', f))
    close_system(sprintf('casper_hw_test_au50_16k_%smhz', f))
    
    disp(sprintf('casper_hw_test_au50_16k_%smhz DONE', f))
    disp('Waiting for 5 minutes')
    pause(300)

    load_system('casper_hw_test_au50_64k')
    set_param('casper_hw_test_au50_64k/AU50', 'clk_rate', f)
    save_system('casper_hw_test_au50_64k', sprintf('./fft_testing/hardware/casper/4in/casper_hw_test_au50_64k_%smhz', f))
    set_param(sprintf('casper_hw_test_au50_64k_%smhz', f),'SimulationCommand','Update')
    jasper(sprintf('casper_hw_test_au50_64k_%smhz', f))
    save_system(sprintf('casper_hw_test_au50_64k_%smhz', f))
    close_system(sprintf('casper_hw_test_au50_64k_%smhz', f))
    
    disp(sprintf('casper_hw_test_au50_64k_%smhz DONE', f))
    disp('Waiting for 5 minutes')
    pause(300)
    
    %******************** ASTRON ********************
    load_system('astron_hw_test_au50_1k')
    set_param('astron_hw_test_au50_1k/AU50', 'clk_rate', f)
    save_system('astron_hw_test_au50_1k', sprintf('./fft_testing/hardware/astron/4in/astron_hw_test_au50_1k_%smhz', f))
    set_param(sprintf('astron_hw_test_au50_1k_%smhz', f),'SimulationCommand','Update')
    jasper(sprintf('astron_hw_test_au50_1k_%smhz', f))
    save_system(sprintf('astron_hw_test_au50_1k_%smhz', f))
    close_system(sprintf('astron_hw_test_au50_1k_%smhz', f))
    
    disp(sprintf('astron_hw_test_au50_1k_%smhz DONE', f))
    disp('Waiting for 5 minutes')
    pause(300)

    load_system('astron_hw_test_au50_16k')
    set_param('astron_hw_test_au50_16k/AU50', 'clk_rate', f)
    save_system('astron_hw_test_au50_16k', sprintf('./fft_testing/hardware/astron/4in/astron_hw_test_au50_16k_%smhz', f))
    set_param(sprintf('astron_hw_test_au50_16k_%smhz', f),'SimulationCommand','Update')
    jasper(sprintf('astron_hw_test_au50_16k_%smhz', f))
    save_system(sprintf('astron_hw_test_au50_16k_%smhz', f))
    close_system(sprintf('astron_hw_test_au50_16k_%smhz', f))
    
    disp(sprintf('astron_hw_test_au50_16k_%smhz DONE', f))
    disp('Waiting for 5 minutes')
    pause(300)

    %******************** XILINX ********************
    load_system('xilinx_ssr_hw_test_au50_1k')
    set_param('xilinx_ssr_hw_test_au50_1k/AU50', 'clk_rate', f)
    save_system('xilinx_ssr_hw_test_au50_1k', sprintf('./fft_testing/hardware/xssr/4in/xilinx_ssr_hw_test_au50_1k_%smhz', f))
    set_param(sprintf('xilinx_ssr_hw_test_au50_1k_%smhz', f),'SimulationCommand','Update')
    jasper(sprintf('xilinx_ssr_hw_test_au50_1k_%smhz', f))
    save_system(sprintf('xilinx_ssr_hw_test_au50_1k_%smhz', f))
    close_system(sprintf('xilinx_ssr_hw_test_au50_1k_%smhz', f))
    
    disp(sprintf('xilinx_ssr_hw_test_au50_1k_%smhz DONE', f))
    disp('Waiting for 5 minutes')
    pause(300)

    load_system('xilinx_ssr_hw_test_au50_16k')
    set_param('xilinx_ssr_hw_test_au50_16k/AU50', 'clk_rate', f)
    save_system('xilinx_ssr_hw_test_au50_16k', sprintf('./fft_testing/hardware/xssr/4in/xilinx_ssr_hw_test_au50_16k_%smhz', f))
    set_param(sprintf('xilinx_ssr_hw_test_au50_16k_%smhz', f),'SimulationCommand','Update')
    jasper(sprintf('xilinx_ssr_hw_test_au50_16k_%smhz', f))
    save_system(sprintf('xilinx_ssr_hw_test_au50_16k_%smhz', f))
    close_system(sprintf('xilinx_ssr_hw_test_au50_16k_%smhz', f))
    
    disp(sprintf('xilinx_ssr_hw_test_au50_16k_%smhz DONE', f))
    disp('Waiting for 5 minutes')
    pause(300)
    
    load_system('xilinx_ssr_hw_test_au50_64k')
    set_param('xilinx_ssr_hw_test_au50_64k/AU50', 'clk_rate', f)
    save_system('xilinx_ssr_hw_test_au50_64k', sprintf('./fft_testing/hardware/xssr/4in/xilinx_ssr_hw_test_au50_64k_%smhz', f))
    set_param(sprintf('xilinx_ssr_hw_test_au50_64k_%smhz', f),'SimulationCommand','Update')
    jasper(sprintf('xilinx_ssr_hw_test_au50_64k_%smhz', f))
    save_system(sprintf('xilinx_ssr_hw_test_au50_64k_%smhz', f))
    close_system(sprintf('xilinx_ssr_hw_test_au50_64k_%smhz', f))
    
    disp(sprintf('xilinx_ssr_hw_test_au50_64k_%smhz DONE', f))
    disp('Waiting for 5 minutes')
    pause(300)
    
    tEnd = toc(tStart)
 
end

clear

%% AU50 (work/repos/mlib_devel): wideband factor = 8

fs = 800e6;             % sampling frequency
NFFT = 2^10;            % FFT Size
binbw = fs/NFFT *1e-6;  % BW in MHz in a single bin
N = NFFT*2;            % length of time vector
t = (1:N)/fs;           % time vector
a_sig = 0.45;         % tone amplitude
f_sig = 50*1e6;         % tone frequency 

real_sig = a_sig*cos(2*pi*f_sig*t);

sim_len = N;
SSR = 4;
bram_depth = N/SSR;
bram_width = log2(bram_depth);

% demux shared_bram
d0 = int32(real_sig(1:8:end)*(2^15));
d1 = int32(real_sig(2:8:end)*(2^15));
d2 = int32(real_sig(3:8:end)*(2^15));
d3 = int32(real_sig(4:8:end)*(2^15));
d4 = int32(real_sig(5:8:end)*(2^15));
d5 = int32(real_sig(6:8:end)*(2^15));
d6 = int32(real_sig(7:8:end)*(2^15));
d7 = int32(real_sig(8:8:end)*(2^15));

freqs = [600, 550, 500, 450, 400, 250, 200];

for i = 1:1
    f = string(freqs(i));

    %******************** CASPER ********************
    load_system('casper_hw_test_au50_1k_8in')
    set_param('casper_hw_test_au50_1k_8in/AU50', 'clk_rate', f)
    save_system('casper_hw_test_au50_1k_8in', sprintf('./fft_testing/hardware/casper/8in/casper_hw_test_au50_1k_8in_%smhz', f))
    set_param(sprintf('casper_hw_test_au50_1k_8in_%smhz', f),'SimulationCommand','Update')
\      save_system()
    jasper(sprintf('casper_hw_test_au50_1k_8in_%smhz', f))
    save_system(sprintf('casper_hw_test_au50_1k_8in_%smhz', f))
    close_system(sprintf('casper_hw_test_au50_1k_8in_%smhz', f))

    load_system('casper_hw_test_au50_16k_8in')
    set_param('casper_hw_test_au50_16k_8in/AU50', 'clk_rate', f)
    save_system('casper_hw_test_au50_16k_8in', sprintf('./fft_testing/hardware/casper/8in/casper_hw_test_au50_16k_8in_%smhz', f))
    set_param(sprintf('casper_hw_test_au50_16k_8in_%smhz', f),'SimulationCommand','Update')
    save_system()
    jasper(sprintf('casper_hw_test_au50_16k_8in_%smhz', f))
    save_system(sprintf('casper_hw_test_au50_16k_8in_%smhz', f))
    close_system(sprintf('casper_hw_test_au50_16k_8in_%smhz', f))

    load_system('casper_hw_test_au50_64k_8in')
    set_param('casper_hw_test_au50_64k_8in/AU50', 'clk_rate', f)
    save_system('casper_hw_test_au50_64k_8in', sprintf('./fft_testing/hardware/casper/8in/casper_hw_test_au50_64k_8in_%smhz', f))
    set_param(sprintf('casper_hw_test_au50_64k_8in_%smhz', f),'SimulationCommand','Update')
    save_system()
    jasper(sprintf('casper_hw_test_au50_64k_8in_%smhz', f))
    save_system(sprintf('casper_hw_test_au50_64k_8in_%smhz', f))
    close_system(sprintf('casper_hw_test_au50_64k_8in_%smhz', f))

    %******************** ASTRON ********************
    load_system('astron_hw_test_au50_1k_8in')
    set_param('astron_hw_test_au50_1k_8in/AU50', 'clk_rate', f)
    save_system('astron_hw_test_au50_1k_8in', sprintf('./fft_testing/hardware/astron/8in/astron_hw_test_au50_1k_8in_%smhz', f))
    set_param(sprintf('astron_hw_test_au50_1k_8in_%smhz', f),'SimulationCommand','Update')
    save_system()
    jasper(sprintf('astron_hw_test_au50_1k_8in_%smhz', f))
    save_system(sprintf('astron_hw_test_au50_1k_8in_%smhz', f))
    close_system(sprintf('astron_hw_test_au50_1k_8in_%smhz', f))

    load_system('astron_hw_test_au50_16k_8in')
    set_param('astron_hw_test_au50_16k_8in/AU50', 'clk_rate', f)
    save_system('astron_hw_test_au50_16k_8in', sprintf('./fft_testing/hardware/astron/8in/astron_hw_test_au50_16k_8in_%smhz', f))
    set_param(sprintf('astron_hw_test_au50_16k_8in_%smhz', f),'SimulationCommand','Update')
    save_system()
    jasper(sprintf('astron_hw_test_au50_16k_8in_%smhz', f))
    save_system(sprintf('astron_hw_test_au50_16k_8in_%smhz', f))
    close_system(sprintf('astron_hw_test_au50_16k_8in_%smhz', f))

    %******************** XILINX ********************
    load_system('xilinx_ssr_hw_test_au50_1k_8in')
    set_param('xilinx_ssr_hw_test_au50_1k_8in/AU50', 'clk_rate', f)
    save_system('xilinx_ssr_hw_test_au50_1k_8in', sprintf('./fft_testing/hardware/xssr/8in/xilinx_ssr_hw_test_au50_1k_8in_%smhz', f))
    set_param(sprintf('xilinx_ssr_hw_test_au50_1k_8in_%smhz', f),'SimulationCommand','Update')
    save_system()
    jasper(sprintf('xilinx_ssr_hw_test_au50_1k_8in_%smhz', f))
    save_system(sprintf('xilinx_ssr_hw_test_au50_1k_8in_%smhz', f))
    close_system(sprintf('xilinx_ssr_hw_test_au50_1k_8in_%smhz', f))

    load_system('xilinx_ssr_hw_test_au50_16k_8in')
    set_param('xilinx_ssr_hw_test_au50_16k_8in/AU50', 'clk_rate', f)
    save_system('xilinx_ssr_hw_test_au50_16k_8in', sprintf('./fft_testing/hardware/xssr/8in/xilinx_ssr_hw_test_au50_16k_8in_%smhz', f))
    set_param(sprintf('xilinx_ssr_hw_test_au50_16k_8in_%smhz', f),'SimulationCommand','Update')
    save_system()
    jasper(sprintf('xilinx_ssr_hw_test_au50_16k_8in_%smhz', f))
    save_system(sprintf('xilinx_ssr_hw_test_au50_16k_8in_%smhz', f))
    close_system(sprintf('xilinx_ssr_hw_test_au50_16k_8in_%smhz', f))
    
    load_system('xilinx_ssr_hw_test_au50_64k_8in')
    set_param('xilinx_ssr_hw_test_au50_64k_8in/AU50', 'clk_rate', f)
    save_system('xilinx_ssr_hw_test_au50_64k_8in', sprintf('./fft_testing/hardware/xssr/8in/xilinx_ssr_hw_test_au50_64k_8in_%smhz', f))
    set_param(sprintf('xilinx_ssr_hw_test_au50_64k_8in_%smhz', f),'SimulationCommand','Update')
    save_system()
    jasper(sprintf('xilinx_ssr_hw_test_au50_64k_8in_%smhz', f))
    save_system(sprintf('xilinx_ssr_hw_test_au50_64k_8in_%smhz', f))
    close_system(sprintf('xilinx_ssr_hw_test_au50_64k_8in_%smhz', f))
 
end

%% SKARAB (university/test_system/mlib_devel): wideband factor = 4

% %**** CASPER ****
% load_system('casper_hw_test_skarab_1k')
% set_param('casper_hw_test_skarab_1k','SimulationCommand','Update')
% jasper('casper_hw_test_skarab_1k')
% save_system()
% close_system('casper_hw_test_skarab_1k')
% 
% load_system('casper_hw_test_skarab_16k')
% set_param('casper_hw_test_skarab_16k','SimulationCommand','Update')
% jasper('casper_hw_test_skarab_16k')
% save_system()
% close_system('casper_hw_test_skarab_16k')
% 
% load_system('casper_hw_test_skarab_64k')
% set_param('casper_hw_test_skarab_64k','SimulationCommand','Update')
% jasper('casper_hw_test_skarab_64k')
% save_system()
% close_system('casper_hw_test_skarab_64k')
% 
% %**** ASTRON ****
% load_system('astron_hw_test_skarab_1k')
% set_param('astron_hw_test_skarab_1k','SimulationCommand','Update')
% jasper('astron_hw_test_skarab_1k')
% save_system()
% close_system('astron_hw_test_skarab_1k')
% 
% load_system('astron_hw_test_skarab_16k')
% set_param('astron_hw_test_skarab_16k','SimulationCommand','Update')
% jasper('astron_hw_test_skarab_16k')
% save_system()
% close_system('astron_hw_test_skarab_16k')

%% SKARAB (university/test_system/mlib_devel): wideband factor = 8

%**** CASPER ****
% load_system('casper_hw_test_skarab_1k_8in')
% set_param('casper_hw_test_skarab_1k_8in','SimulationCommand','Update')
% jasper('casper_hw_test_skarab_1k_8in')
% save_system()
% close_system('casper_hw_test_skarab_1k_8in')
% 
% load_system('casper_hw_test_skarab_16k_8in')
% set_param('casper_hw_test_skarab_16k_8in','SimulationCommand','Update')
% jasper('casper_hw_test_skarab_16k_8in')
% save_system()
% close_system('casper_hw_test_skarab_16k_8in')
% 
% load_system('casper_hw_test_skarab_64k_8in')
% set_param('casper_hw_test_skarab_64k_8in','SimulationCommand','Update')
% jasper('casper_hw_test_skarab_64k_8in')
% save_system()
% close_system('casper_hw_test_skarab_64k_8in')
% 
% %**** ASTRON ****
% load_system('astron_hw_test_skarab_1k_8in')
% set_param('astron_hw_test_skarab_1k_8in','SimulationCommand','Update')
% jasper('astron_hw_test_skarab_1k_8in')
% save_system()
% close_system('astron_hw_test_skarab_1k_8in')
% 
% load_system('astron_hw_test_skarab_16k_8in')
% set_param('astron_hw_test_skarab_16k_8in','SimulationCommand','Update')
% jasper('astron_hw_test_skarab_16k_8in')
% save_system()
% close_system('astron_hw_test_skarab_16k_8in')



