%% AU50 (work/repos/mlib_devel): wideband factor = 4

% fs = 800e6;             % sampling frequency
% NFFT = 2^16;            % FFT Size
% N = NFFT*2;            % length of time vector
% t = (1:N)/fs;           % time vector
% a_sig = 0.45;         % tone amplitude
% f_sig = 50*1e6;         % tone frequency 
% 
% real_sig = a_sig*cos(2*pi*f_sig*t);
% 
% sim_len = N;
% SSR = 4;
% bram_depth = N/SSR;
% bram_width = log2(bram_depth);
% sr = 1023;
% 
% % demux shared_bram
% d0 = int32(real_sig(1:4:end)*(2^15));
% d1 = int32(real_sig(2:4:end)*(2^15));
% d2 = int32(real_sig(3:4:end)*(2^15));
% d3 = int32(real_sig(4:4:end)*(2^15));
% 
% freqs = [160];
% 
% for i = 1:size(freqs,2)
%     f = string(freqs(i));
% 
%     addpath('/media/morag/linux_storage/storage_home/university/test_system/fft_testing/hardware/builds/designs/4in/au50/templates')
%     filelist = dir('/media/morag/linux_storage/storage_home/university/test_system/fft_testing/hardware/builds/designs/4in/au50/templates')
%     filelist = filelist(~startsWith({filelist.name}, 'tmp'));
%     template_designs = strings;
% 
%     for p = 3:(size(filelist,1))
%         template_designs(p-2) = filelist(p).name;
%     end
% 
%     logfile = sprintf('/media/morag/linux_storage/storage_home/university/test_system/fft_testing/hardware/builds/logs/au50_4in_%smhz_%s.log', f, datestr(now,'yyyy_mm_dd_HH_MM'))
%     fileID = fopen(logfile, 'w')
%     fprintf(fileID, sprintf('%s: Starting %sMHz builds for AU50 (4 inputs)\n\n', datestr(now,'yyyy:mm:dd HH:MM'), f))
% 
%     tStart = tic;
%     for j = 1:(size(template_designs,2))
%         template_design = erase(template_designs(j), '.slx');
%         au50_build_4in(fileID, template_design, f)
% 
%     end
% 
%     fprintf(fileID, '---------------------------------------------------------------\n\n');
% 
%     tEnd = toc(tStart);
%     fprintf(fileID, sprintf('%s: Finished %sMHz builds for AU50 (4 inputs)\n\n', datestr(now,'yyyy:mm:dd HH:MM'), f));
%     fprintf(fileID, sprintf('Total time elapsed: %g hours %g minutes, %g seconds\nn', floor(tEnd/3600), floor(mod(tEnd, 3600)/60), mod(tEnd,60)));
% 
%     fprintf(fileID, '---------------------------------------------------------------\n\n');
%     fclose(fileID);
% end

%% AU50 (work/repos/mlib_devel): wideband factor = 8
clear

fs = 800e6;             % sampling frequency
NFFT = 2^16;            % FFT Size
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
sr = 1023;

% demux shared_bram
d0 = int32(real_sig(1:8:end)*(2^15));
d1 = int32(real_sig(2:8:end)*(2^15));
d2 = int32(real_sig(3:8:end)*(2^15));
d3 = int32(real_sig(4:8:end)*(2^15));
d4 = int32(real_sig(5:8:end)*(2^15));
d5 = int32(real_sig(6:8:end)*(2^15));
d6 = int32(real_sig(7:8:end)*(2^15));
d7 = int32(real_sig(8:8:end)*(2^15));

freqs = [160];

for i = 1:size(freqs,2)
    f = string(freqs(i));

    addpath('/media/morag/linux_storage/storage_home/university/test_system/fft_testing/hardware/builds/designs/8in/au50/templates')
    filelist = dir('/media/morag/linux_storage/storage_home/university/test_system/fft_testing/hardware/builds/designs/8in/au50/templates')
    filelist = filelist(~startsWith({filelist.name}, 'tmp'));
    template_designs = strings;

    for p = 3:(size(filelist,1))
        template_designs(p-2) = filelist(p).name;
    end

    logfile = sprintf('/media/morag/linux_storage/storage_home/university/test_system/fft_testing/hardware/builds/logs/au50_8in_%smhz_%s.log', f, datestr(now,'yyyy_mm_dd_HH_MM'))
    fileID = fopen(logfile, 'w')
    fprintf(fileID, sprintf('%s: Starting %sMHz builds for AU50 (8 inputs)\n\n', datestr(now,'yyyy:mm:dd HH:MM'), f));

    tStart = tic;
    for j = 1:(size(template_designs,2))
        template_design = erase(template_designs(j), '.slx');
        au50_build_8in(fileID, template_design, f)

    end

    fprintf(fileID, '---------------------------------------------------------------\n\n');

    tEnd = toc(tStart);
    fprintf(fileID, sprintf('%s: Finished %sMHz builds for AU50 (8 inputs)\n\n', datestr(now,'yyyy:mm:dd HH:MM'), f));
    fprintf(fileID, sprintf('Total time elapsed: %g hours %g minutes, %g seconds\nn', floor(tEnd/3600), floor(mod(tEnd, 3600)/60), mod(tEnd,60)));

    fprintf(fileID, '---------------------------------------------------------------\n\n');
    fclose(fileID);
end

%% SKARAB (university/test_system/mlib_devel): wideband factor = 4

% clear
% 
% fs = 800e6;             % sampling frequency
% NFFT = 2^10;            % FFT Size
% N = NFFT*2;            % length of time vector
% t = (1:N)/fs;           % time vector
% a_sig = 0.45;         % tone amplitude
% f_sig = 50*1e6;         % tone frequency 
% 
% real_sig = a_sig*cos(2*pi*f_sig*t);
% 
% sim_len = N;
% SSR = 4;
% bram_depth = N/SSR;
% bram_width = log2(bram_depth);
% sr = 1023;
% 
% % demux shared_bram
% d0 = int32(real_sig(1:4:end)*(2^15));
% d1 = int32(real_sig(2:4:end)*(2^15));
% d2 = int32(real_sig(3:4:end)*(2^15));
% d3 = int32(real_sig(4:4:end)*(2^15));
% 
% freqs = [160];
% 
% for i = 1:size(freqs,2)
%     f = string(freqs(i));
% 
%     addpath('/media/morag/linux_storage/storage_home/university/test_system/fft_testing/hardware/builds/designs/4in/skarab/templates')
%     filelist = dir('/media/morag/linux_storage/storage_home/university/test_system/fft_testing/hardware/builds/designs/4in/skarab/templates')
%     filelist = filelist(~startsWith({filelist.name}, 'tmp'));
%     template_designs = strings;
% 
%     for p = 3:(size(filelist,1))
%         template_designs(p-2) = filelist(p).name;
%     end
% 
%     logfile = sprintf('/media/morag/linux_storage/storage_home/university/test_system/fft_testing/hardware/builds/logs/skarab_4in_%smhz_%s.log', f, datestr(now,'yyyy_mm_dd_HH_MM'))
%     fileID = fopen(logfile, 'w')
%     fprintf(fileID, sprintf('%s: Starting %sMHz builds for SKARAB (4 inputs)\n\n', datestr(now,'yyyy:mm:dd HH:MM'), f))
% 
%     tStart = tic;
%     for j = 1:(size(template_designs,2))
%         template_design = erase(template_designs(j), '.slx');
%         skarab_build_4in(fileID, template_design, f)
% 
%     end
% 
%     fprintf(fileID, '---------------------------------------------------------------\n\n');
% 
%     tEnd = toc(tStart);
%     fprintf(fileID, sprintf('%s: Finished %sMHz builds for SKARAB (4 inputs)\n\n', datestr(now,'yyyy:mm:dd HH:MM'), f));
%     fprintf(fileID, sprintf('Total time elapsed: %g hours %g minutes, %g seconds\n\n', floor(tEnd/3600), floor(mod(tEnd, 3600)/60), mod(tEnd,60)));
% 
%     fprintf(fileID, '---------------------------------------------------------------\n\n');
%     fclose(fileID);
% end

%% SKARAB (university/test_system/mlib_devel): wideband factor = 8

% clear
% 
% fs = 800e6;             % sampling frequency
% NFFT = 2^10;            % FFT Size
% binbw = fs/NFFT *1e-6;  % BW in MHz in a single bin
% N = NFFT*2;            % length of time vector
% t = (1:N)/fs;           % time vector
% a_sig = 0.45;         % tone amplitude
% f_sig = 50*1e6;         % tone frequency 
% 
% real_sig = a_sig*cos(2*pi*f_sig*t);
% 
% sim_len = N;
% SSR = 4;
% bram_depth = N/SSR;
% bram_width = log2(bram_depth);
% sr = 1023;
% 
% % demux shared_bram
% d0 = int32(real_sig(1:8:end)*(2^15));
% d1 = int32(real_sig(2:8:end)*(2^15));
% d2 = int32(real_sig(3:8:end)*(2^15));
% d3 = int32(real_sig(4:8:end)*(2^15));
% d4 = int32(real_sig(5:8:end)*(2^15));
% d5 = int32(real_sig(6:8:end)*(2^15));
% d6 = int32(real_sig(7:8:end)*(2^15));
% d7 = int32(real_sig(8:8:end)*(2^15));
% 
% freqs = [160];
% 
% for i = 1:size(freqs,2)
%     f = string(freqs(i));
% 
%     addpath('/media/morag/linux_storage/storage_home/university/test_system/fft_testing/hardware/builds/designs/8in/skarab/templates')
%     filelist = dir('/media/morag/linux_storage/storage_home/university/test_system/fft_testing/hardware/builds/designs/8in/skarab/templates')
%     filelist = filelist(~startsWith({filelist.name}, 'tmp'));
%     template_designs = strings;
% 
%     for p = 3:(size(filelist,1))
%         template_designs(p-2) = filelist(p).name;
%     end
% 
%     logfile = sprintf('/media/morag/linux_storage/storage_home/university/test_system/fft_testing/hardware/builds/logs/skarab_8in_%smhz_%s.log', f, datestr(now,'yyyy_mm_dd_HH_MM'))
%     fileID = fopen(logfile, 'w')
%     fprintf(fileID, sprintf('%s: Starting %sMHz builds for SKARAB (8 inputs)\n\n', datestr(now,'yyyy:mm:dd HH:MM'), f))
% 
%     tStart = tic;
%     for j = 1:(size(template_designs,2))
%         template_design = erase(template_designs(j), '.slx');
%         skarab_build_8in(fileID, template_design, f)
% 
%     end
% 
%     fprintf(fileID, '---------------------------------------------------------------\n\n');
% 
%     tEnd = toc(tStart);
%     fprintf(fileID, sprintf('%s: Finished %sMHz builds for SKARAB (8 inputs)\n\n', datestr(now,'yyyy:mm:dd HH:MM'), f));
%     fprintf(fileID, sprintf('Total time elapsed: %g hours %g minutes, %g seconds\nn', floor(tEnd/3600), floor(mod(tEnd, 3600)/60), mod(tEnd,60)));
% 
%     fprintf(fileID, '---------------------------------------------------------------\n\n')
%     fclose(fileID);
% end




