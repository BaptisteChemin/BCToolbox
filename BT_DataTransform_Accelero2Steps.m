

clear;
% load the .mat file with accelero signals
[matnames, matpath] = uigetfile('*.mat', 'Get the stimulation .wav file','MultiSelect','on');
if ischar (matnames);matnames = {char(matnames)}; end



clc; 
for f=1:length(matnames)
    matname = char(matnames(f));
    disp (matname);
    eval(sprintf('%s',['load(''' char(matpath) char(matname) ''');']));
    
    
    
    
    trigger = trigger/max(abs(trigger));
    
    fs=250;
    
    
    % get the best channel for both legs, and get trigger channel
    
    
    figure; set(gcf,'units','normalized','outerposition',[0 0 1 1])
    hold on
    plot(jambe_gauche_x,'color',[.8 .8 .8]);
    plot(trigger*1000-500,'color',[.95 .95 .95]);
    plot(jambe_droite_x,'b');
    plot(jambe_droite_y,'r');
    legend('jambe gauche','triggers','jambe droite x', 'jambe droite y')
    hold off
    xlim([length(jambe_gauche_x)/4 2*length(jambe_gauche_x)/4]);
    ylim([min([jambe_droite_x jambe_droite_y]) max([jambe_droite_x jambe_droite_y])])
    
    choice = questdlg('Which is the best channel?','Channel Selection','x','y','x');
    eval(sprintf('%s',['jambe_droite = jambe_droite_' char(choice) ';']));
    close
    
    figure; set(gcf,'units','normalized','outerposition',[0 0 1 1])
    hold on
    plot(jambe_droite,'color',[.8 .8 .8]);
    plot(trigger*1000-500,'color',[.95 .95 .95]);
    plot(jambe_gauche_x,'b');
    plot(jambe_gauche_y,'r');
    legend('jambe droite','triggers','jambe gauche x', 'jambe gauche y')
    hold off
    xlim([length(jambe_gauche_x)/4 2*length(jambe_gauche_x)/4]);
    ylim([min([jambe_gauche_x jambe_gauche_y]) max([jambe_gauche_x jambe_gauche_y])])
    
    
    choice = questdlg('Which is the best channel?','Channel Selection','x','y','x');
    eval(sprintf('%s',['jambe_gauche = jambe_gauche_' char(choice) ';']));
    close
    
    clear jambe_droite_x; clear jambe_droite_y; clear jambe_gauche_x; clear jambe_gauche_y;
    
    plot(jambe_droite); set(gcf,'units','normalized','outerposition',[0 0 1 1])
    xlim([length(jambe_droite)/4 2*length(jambe_droite)/4]);
    direction = questdlg('what is the direction of maximal acceleration?','Channel Orientation','falling','rising','falling');
    if strcmp(direction,'rising')
        jambe_droite = -jambe_droite;
    end
    close
    
    plot(jambe_gauche); set(gcf,'units','normalized','outerposition',[0 0 1 1])
    xlim([length(jambe_droite)/4 2*length(jambe_droite)/4]);
    direction = questdlg('what is the direction of maximal acceleration?','Channel Orientation','falling','rising','falling');
    if strcmp(direction,'rising')
        jambe_gauche = -jambe_gauche;
    end
    close
    
    
    % trigger latencies
    threshold = .25;
    idxl = trigger>= threshold;
    [~,Time_Trigs] = risetime (double(idxl),fs);
    
    
    % identify the max negative acceleration (=step latencies)
    line                = jambe_droite;
    line                = double(line);
    line                = line-nanmean(line);
    line                = line/nanstd(line);
    line(isnan(line))   = nanmean(line);
    
    line_oppo           = jambe_gauche;
    line_oppo           = double(line_oppo);
    line_oppo           = line_oppo-nanmean(line_oppo);
    line_oppo           = line_oppo/nanstd(line_oppo);
    
    
    time = 1/fs:1/fs:length(line)/fs;
    
    
    BT_PreProcess_GetStepLatencies
    Time_jambe_droite = steps_t;
    clear line; clear line_oppo;
    
    line                = jambe_gauche;
    line                = double(line);
    line                = line-nanmean(line);
    line                = line/nanstd(line);
    line(isnan(line))   = nanmean(line);
    
    line_oppo           = jambe_droite;
    line_oppo           = double(line_oppo);
    line_oppo           = line_oppo-nanmean(line_oppo);
    line_oppo           = line_oppo/nanstd(line_oppo);
    
    BT_PreProcess_GetStepLatencies
    Time_jambe_gauche = steps_t;
    
    
    %% save
    eval(sprintf('%s',[ char(matname(1:end-4)) '.trigger_times = Time_Trigs;']));
    eval(sprintf('%s',[ char(matname(1:end-4)) '.jambe_droite_times = Time_jambe_droite;']));
    eval(sprintf('%s',[ char(matname(1:end-4)) '.jambe_gauche_times = Time_jambe_gauche;']));
    
    
    eval(sprintf('%s',['save ([''' char(matname(1:end-4)) '_latencies.mat''], ''' char(matname(1:end-4)) ''')']));
    disp('*** finished ***')
    
end

disp('*** all finished ***')