%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sim_ffts.m
% Author: Morag Brown
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inputs_dir = "./inputs";
designs_dir = "./designs/8in";
outputs_dir = "./outputs/8in";

[input_files, design_files, demux] = get_input_files(inputs_dir, designs_dir);

% Begin simulation for each input file
for i = 1:size(input_files, 2)
    input_file = char(input_files(i));
    fprintf("Using Test Vector: %s\n", input_file);

    % loop through each design
    for k = 1:size(design_files, 2)
        design_file = char(design_files(k));
        fprintf("Loading design: %s\n", design_file);
        load_system(erase(design_file, '.slx'))
        [fft_type, hw_type, nfft] = parse_design_file(design_file);
        
        % loop through each shifting schedule - shift values change according to fft length:
        % shifts are: [none, every stage for last half, every stage for first half, alt stages, all stages]
    %     shift_type = ["none", "lastH", "firstH", "alt", "all"];
    %     shift_dict = containers.Map({2^10, 2^14, 2^16}, ...
    %     {[0, 31, 992, 682, 1023], [0, 127, 16256, 10922, 16383], [0, 255, 65280, 43690, 65535]});
tha
        shift_type = ["none"];
        shift_dict = containers.Map({2^10, 2^14, 2^16}, ...
        {[0], [0], [0]});
    
        shift_schedules = shift_dict(nfft); 
        for j = 1:size(shift_schedules, 2)
            sr_type = shift_type(j);
            sr = shift_schedules(j);
            sim_len = 2*nfft;
            
            [sig, fs, output_file] = parse_input_file(input_file, outputs_dir, fft_type, hw_type, nfft, sr, sr_type, demux);

            % equivalent of 'write_signal' function
            if demux == 4
                d0 = int32(sig(1:4:end)*(2^15)).';
                d1 = int32(sig(2:4:end)*(2^15)).';
                d2 = int32(sig(3:4:end)*(2^15)).';
                d3 = int32(sig(4:4:end)*(2^15)).';

            elseif demux == 8
                d0 = int32(sig(1:8:end)*(2^15)).';
                d1 = int32(sig(2:8:end)*(2^15)).';
                d2 = int32(sig(3:8:end)*(2^15)).';
                d3 = int32(sig(4:8:end)*(2^15)).';
                d4 = int32(sig(5:8:end)*(2^15)).';
                d5 = int32(sig(6:8:end)*(2^15)).';
                d6 = int32(sig(7:8:end)*(2^15)).';
                d7 = int32(sig(8:8:end)*(2^15)).';
            end
            
            fprintf("Saving and Updating Model\n");
            save_system(erase(design_file, '.slx'))
            set_param(erase(design_file, '.slx'),'SimulationCommand','Update')

            fprintf("Beginning Simulation using:\nTest Vector - %s\nDesign File - %s\nShift Type - %s\n", ...
                    input_file, design_file, sr_type);
            tic;
            f_out = sim(erase(design_file, '.slx'), sim_len);
            T = toc;

            output_fft = get_output_data(T, sim_len, demux, fft_type, nfft, output_file, fs, f_out);
            
            % clean up directory and workspace
            delete *.log
            delete *.slxc
            
            clear sr output_file f_out val_id T  
            
        end       
        save_system(erase(design_file, '.slx'))
        close_system(erase(design_file, '.slx'))   
        clear design_file 
    end
    clear input_file 
end

function [input_files, design_files, demux] = get_input_files(inputs_dir, designs_dir)
    inputs_list = dir(inputs_dir);
    designs_list = dir(designs_dir);

    input_files = {inputs_list(~[inputs_list.isdir]).name};
    design_files = {designs_list(~[designs_list.isdir]).name};
    
    demux_split = split(designs_dir, '/');
    demux_dict = containers.Map({'4in', '8in'}, [4, 8]);
    demux = demux_dict(demux_split(end));
     
end

function [fft_type, hw_type, nfft] = parse_design_file(design_file)
    split_design = split(design_file, '_');
    fft_type = char(split_design(1));
    hw_type = char(split_design(4));
    length_split = split(split_design(5), '.');
    fft_length = char(length_split(1));
    
    length_dict = containers.Map({'1k', '16k', '64k'}, [2^10, 2^14, 2^16]);
    nfft = length_dict(fft_length);
end

function [sig, fs, output_file] = parse_input_file(input_file, outputs_dir, fft_type, hw_type, nfft, sr, sr_type, demux)
    file = fopen(input_file);
    test_vector = fscanf(file, '%f');
    sig = test_vector(1:nfft);
    
    fs_split = split(input_file, '_');
    i = find(contains(fs_split,'msps'));
    pat = digitsPattern;
    fs_str = char(extract(fs_split(i),pat));
    fs = str2num(fs_str);
    
    output_file = outputs_dir + "/" + hw_type + "_" + fft_type + "_" + num2str(nfft) + "_" + input_file + "_sr-" + sr_type + ".output";
end

function output_fft = get_output_data(T, sim_len, demux, fft_type, nfft, output_file, fs, f_out)

    if fft_type == 'casper'
        if demux == 4
            dmux_out = 2;

            % find valid data index, sync_out is a pulse
            % output data valid 1 clock cycle after sync_out goes high
            % add 1 to sim len cos matlab indexes from 1, not 0 (eg 1:2049, not 0:2048)
            % add dmux_out to val_len also because of indexing from 1 (+1 to each
            % val_len per demuxed lane)
            val_id = find(f_out.sync_out) + 1;
            val_len = (sim_len + 1)*dmux_out - val_id*dmux_out + dmux_out;

            % interleave output
            fft_re(1:2:val_len) = f_out.out_re(val_id:end);
            fft_re(2:2:val_len) = f_out.out_re1(val_id:end);

            fft_im(1:2:val_len) = f_out.out_im(val_id:end);
            fft_im(2:2:val_len) = f_out.out_im1(val_id:end);

            fft_complex = fft_re(1:nfft/2) + fft_im(1:nfft/2)*1j;  % CASPER already discards neg freq components, (N/2 : N) is next spectra
            output_fft = abs(fft_complex);

        elseif demux == 8
            dmux_out = 4;

            % find valid data index, sync_out is a pulse
            % output data valid 1 clock cycle after sync_out goes high
            % add 1 to sim len cos matlab indexes from 1, not 0 (eg 1:2049, not 0:2048)
            % add dmux_out to val_len also because of indexing from 1 (+1 to each
            % val_len per demuxed lane)
            val_id = find(f_out.sync_out) + 1;
            val_len = (sim_len + 1)*dmux_out - val_id*dmux_out + dmux_out;

            % interleave output
            fft_re(1:4:val_len) = f_out.out_re(val_id:end);
            fft_re(2:4:val_len) = f_out.out_re1(val_id:end);
            fft_re(3:4:val_len) = f_out.out_re2(val_id:end);
            fft_re(4:4:val_len) = f_out.out_re3(val_id:end);

            fft_im(1:4:val_len) = f_out.out_im(val_id:end);
            fft_im(2:4:val_len) = f_out.out_im1(val_id:end);
            fft_im(3:4:val_len) = f_out.out_im2(val_id:end);
            fft_im(4:4:val_len) = f_out.out_im3(val_id:end);

            fft_complex = fft_re(1:nfft/2) + fft_im(1:nfft/2)*1j;  % CASPER already discards neg freq components, (N/2 : N) is next spectra
            output_fft = abs(fft_complex);
        end    
        
    
    elseif fft_type == 'astron'
        if demux == 4
            dmux_out = 4;

            % find valid data index, sync_out is a pulse
            % output data valid 1 clock cycle after sync_out goes high
            % add 1 to sim len cos matlab indexes from 1, not 0 (eg 1:2049, not 0:2048)
            % add dmux_out to val_len also because of indexing from 1 (+1 to each
            % val_len per demuxed lane)
            % divide val_len by 2 cos spectrum of 2 signals in each output lane
            val_id = find(f_out.sync_out) + 1;
            val_len = ((sim_len + 1)*dmux_out - val_id*dmux_out + dmux_out)/2;
            
            % interleave output
            fft_re(1:nfft/8) = f_out.out_re(val_id:2:val_id + nfft/4 - 1);
            fft_re(nfft/8 + 1:2*nfft/8) = f_out.out_re1(val_id:2:val_id + nfft/4 - 1);
            fft_re(2*nfft/8 + 1:3*nfft/8) = f_out.out_re2(val_id:2:val_id + nfft/4 - 1);
            fft_re(3*nfft/8 + 1:4*nfft/8) = f_out.out_re3(val_id:2:val_id + nfft/4 - 1);
            
            fft_im(1:nfft/8) = f_out.out_im(val_id:2:val_id + nfft/4 - 1);
            fft_im(nfft/8 + 1:2*nfft/8) = f_out.out_im1(val_id:2:val_id + nfft/4 - 1);
            fft_im(2*nfft/8 + 1:3*nfft/8) = f_out.out_im2(val_id:2:val_id + nfft/4 - 1);
            fft_im(3*nfft/8 + 1:4*nfft/8) = f_out.out_im3(val_id:2:val_id + nfft/4 - 1);
            
            fft_complex = fft_re + fft_im*1j;
            output_fft = abs(fft_complex);
            
        elseif demux == 8
            dmux_out = 8;

            % find valid data index, sync_out is a pulse
            % output data valid 1 clock cycle after sync_out goes high
            % add 1 to sim len cos matlab indexes from 1, not 0 (eg 1:2049, not 0:2048)
            % add dmux_out to val_len also because of indexing from 1 (+1 to each
            % val_len per demuxed lane)
            % divide val_len by 2 cos spectrum of 2 signals in each output lane
            val_id = find(f_out.sync_out) + 1;
            val_len = ((sim_len + 1)*dmux_out - val_id*dmux_out + dmux_out)/2;
            
            % interleave output
            fft_re(1:nfft/16) = f_out.out_re(val_id:2:val_id + nfft/8 - 1);
            fft_re(nfft/16 + 1:2*nfft/16) = f_out.out_re1(val_id:2:val_id + nfft/8 - 1);
            fft_re(2*nfft/16 + 1:3*nfft/16) = f_out.out_re2(val_id:2:val_id + nfft/8 - 1);
            fft_re(3*nfft/16 + 1:4*nfft/16) = f_out.out_re3(val_id:2:val_id + nfft/8 - 1);
            fft_re(4*nfft/16 + 1:5*nfft/16) = f_out.out_re4(val_id:2:val_id + nfft/8 - 1);
            fft_re(5*nfft/16 + 1:6*nfft/16) = f_out.out_re5(val_id:2:val_id + nfft/8 - 1);
            fft_re(6*nfft/16 + 1:7*nfft/16) = f_out.out_re6(val_id:2:val_id + nfft/8 - 1);
            fft_re(7*nfft/16 + 1:8*nfft/16) = f_out.out_re7(val_id:2:val_id + nfft/8 - 1);
            
            fft_im(1:nfft/16) = f_out.out_im(val_id:2:val_id + nfft/8 - 1);
            fft_im(nfft/16 + 1:2*nfft/16) = f_out.out_im1(val_id:2:val_id + nfft/8 - 1);
            fft_im(2*nfft/16 + 1:3*nfft/16) = f_out.out_im2(val_id:2:val_id + nfft/8 - 1);
            fft_im(3*nfft/16 + 1:4*nfft/16) = f_out.out_im3(val_id:2:val_id + nfft/8 - 1);
            fft_im(4*nfft/16 + 1:5*nfft/16) = f_out.out_im4(val_id:2:val_id + nfft/8 - 1);
            fft_im(5*nfft/16 + 1:6*nfft/16) = f_out.out_im5(val_id:2:val_id + nfft/8 - 1);
            fft_im(6*nfft/16 + 1:7*nfft/16) = f_out.out_im6(val_id:2:val_id + nfft/8 - 1);
            fft_im(7*nfft/16 + 1:8*nfft/16) = f_out.out_im7(val_id:2:val_id + nfft/8 - 1);
            
            fft_complex = fft_re + fft_im*1j;
            output_fft = abs(fft_complex);
            
        end 
        
    elseif fft_type == 'xilinx'
        val_id = find(f_out.val_out~=0, 1, 'first');
        % extract rows from val_id to last
        fft_re = f_out.out_re(val_id:end,:);
        fft_im = f_out.out_im(val_id:end,:);

        fft_complex =  fft_re + fft_im*1j;
        fft_complex_flat = reshape((fft_complex).',1,(size(fft_complex,1)*size(fft_complex,2)));
        x_output_fft = abs(fftshift(fft_complex_flat(1:nfft)));

        output_fft = abs(fft_complex_flat(1:nfft/2));
          
    end
    
    if isempty(find(f_out.overflow(val_id:end)))
        of = 0;
    else
        of = 1;
    end
    
    file = fopen(output_file, 'w');
    fprintf(file, 'sim time: %f\n', T);
    fprintf(file, 'latency: %f\n', val_id);
    fprintf(file, 'overflow: %f\n', of);
    fprintf(file, '%e\n', output_fft);
    fclose(file);
    
    clear of
end

