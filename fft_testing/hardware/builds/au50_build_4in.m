function au50_build_4in(fileID, design, clk)
    new_path = sprintf('/media/morag/linux_storage/storage_home/university/test_system/fft_testing/hardware/builds/designs/4in/au50/generated/%s_%smhz', design, clk);
    new_design = sprintf('%s_%smhz', design, clk);
    
    fprintf(fileID, '---------------------------------------------------------------\n\n');
    
    fprintf(fileID, sprintf('%s: Creating %s\n', datestr(now,'yyyy:mm:dd HH:MM'), new_design))    
    load_system(design)
    set_param(sprintf('%s/AU50', design), 'clk_rate', clk)
    save_system(design, new_path)
    
    fprintf(fileID, sprintf('%s: Updating %s\n', datestr(now,'yyyy:mm:dd HH:MM'), new_design))
    set_param(new_design,'SimulationCommand','Update')
    
    fprintf(fileID, sprintf('%s: Building %s\n', datestr(now,'yyyy:mm:dd HH:MM'), new_design))    
    tic
    jasper(new_design)
    save_system(new_design)
    close_system(new_design)     
    T =  toc;
    s = seconds(T);
    s.Format = 'hh:mm:ss';
    
    fprintf(fileID, sprintf('%s: Finished building %s\n', datestr(now,'yyyy:mm:dd HH:MM'), new_design));
    fprintf(fileID, sprintf('Build time: %s\n\n', s));
    
    fprintf(fileID, '---------------------------------------------------------------\n\n');
    
    % wait 5 mins for next build to let model composer release resources
    % (hopefully prevent crashing)
    for k = 5:-1:1
        fprintf('Waiting for %s minutes\n', string(k))
        pause(60)
    end
     
end
    