function [r_wave,lagtime_wave,rMax_wave_Time,rMax_wave,r0_wave,rac_wave,lagactime_wave,r0_Repp,r1_Repp,ac1_Repp] = BT_Analysis_CrossCor(Time_Steps,Time_Trigs,fs,deviants,win,fig)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    %% warning
    warning('Needs Matlab R2017 or higher')
    %% data
    All_Events              = sort([Time_Trigs Time_Steps]);
    time                    = (All_Events(1)-1):1/fs:(All_Events(end)+1);
    
    IOI_Trigs               = diff(Time_Trigs);
    IOI_Steps               = diff(Time_Steps);
    
    [Idx_Steps,Idx_Trigs]   = BT_DataTransform_Time2Idx(Time_Steps,Time_Trigs,time);
    
    %% remove the segments with missing events
    [time_keep,IOI_Trigs_nodeviant,IOI_Steps_nodeviant] = BT_DataTransform_Time2IEI(Time_Steps,Time_Trigs,time,'removedeviants');
    
    %% Adjust the baseline
    IOI_Trigs_bl            = IOI_Trigs-mean(IOI_Trigs_nodeviant);
    IOI_Steps_bl            = IOI_Steps-mean(IOI_Steps_nodeviant);
    
%     THIS WAS USED when the steps of only one leg were used. Now, removing the mean InterEventInterval in each signal should do the job...  
%     if doubletrigs == 1 
%         IOI_Trigs = IOI_Trigs*2;
%     end

    %% limits of analysis
    Idx_Start_Steps_tf      = 2; % Index of the second actual step for time-frequency waveforms
    Idx_Start_Trigs_tf      = find(abs(Time_Trigs-Time_Steps(Idx_Start_Steps_tf))==min(abs(Time_Trigs-Time_Steps(Idx_Start_Steps_tf)))); % Index of the closest event to second step   
    Idx_End_Trigs_tf        = length(Time_Trigs); % Index of the last stimulus
    Idx_End_Steps_tf        = find(abs(Time_Steps-Time_Trigs(Idx_End_Trigs_tf))==min(abs(Time_Steps-Time_Trigs(Idx_End_Trigs_tf)))); % Index of the closest steps to last stimulus
        
    %Time_Start_Steps        = Time_Steps(Idx_Start_Steps_tf); % latency of the second actual step
    %Time_Start_Trigs        = Time_Trigs(Idx_Start_Steps_tf); % latency of the second closest event to second step
    %Time_End_Trigs          = Time_Trigs(Idx_End_Trigs_tf); % latency of the last stimulus
    %Time_End_Steps          = Time_Steps(Idx_End_Steps_tf); % latency of the closest steps to last stimulus
    
    % NOT GOOD! SHOULD DO IT with IOI_x_nodeviant (see above as a start)
    Idx_Start_Steps_IOI     = Idx_Start_Steps_tf - 1; if Idx_Start_Steps_IOI==0; Idx_Start_Steps_IOI=1; end
    Idx_Start_Trigs_IOI     = Idx_Start_Trigs_tf - 1; if Idx_Start_Trigs_IOI==0; Idx_Start_Trigs_IOI=1; end
    Idx_End_Steps_IOI       = Idx_End_Steps_tf - 1;
    Idx_End_Trigs_IOI       = Idx_End_Trigs_tf - 1;
    
    %% cross-correlations on time-IOI waveforms
    % generate the waveforms
    Dots_Trigs = NaN(1,length(time));
    for f=1:length(IOI_Trigs_bl)
        Dots_Trigs(Idx_Trigs(f+1))=IOI_Trigs_bl(f); 
    end    
    Wave_Trigs = fillmissing(Dots_Trigs,'linear');
    
    Dots_Steps = NaN(1,length(time));
    for f=1:length(IOI_Steps_bl)
        Dots_Steps(Idx_Steps(f+1))=IOI_Steps_bl(f);
    end    
    Wave_Steps = fillmissing(Dots_Steps,'linear');
    % remove the missing data
    Wave_Steps_keep = Wave_Steps.*time_keep;
    Wave_Trigs_keep = Wave_Trigs.*time_keep;
    % crop the waveforms so there is no extrapolation due to fillmissing function    
    start       = Idx_Steps(Idx_Start_Steps_tf)+1;
    stop        = Idx_Steps(Idx_End_Steps_tf)-1;
    Wave_Trigs_crop = Wave_Trigs_keep(start:stop);
    Wave_Steps_crop = Wave_Steps_keep(start:stop); 
    
    
        
    Wave_Trigs_ok = Wave_Trigs_crop;
    Wave_Trigs_ok(isnan(Wave_Trigs_ok))=[];
    Wave_Steps_ok = Wave_Steps_crop;
    Wave_Steps_ok(isnan(Wave_Steps_ok))=[];
    % time_ok = time(start:stop) ;
    % time_ok(isnan(time_keep))=[];
    
    [r_wave,lag_wave]=xcorr(Wave_Steps_ok,Wave_Trigs_ok,'coeff');
    lagtime_wave = lag_wave/fs;
    
    window = [find(lagtime_wave==-win) find(lagtime_wave==+win)];
    rMax_wave_Time = lagtime_wave(find(r_wave(window(1):window(2))==max(r_wave(window(1):window(2))))+ window(1));
    rMax_wave = max(r_wave(window(1):window(2)));
    corrmax = lagtime_wave==rMax_wave_Time;
    
    r0_idx = lagtime_wave==0;
    r0_wave = r_wave(r0_idx);
    
    [rac_wave,lagac_wave]=xcorr(Wave_Trigs_ok,Wave_Trigs_ok,'coeff');
    lagactime_wave = lagac_wave/fs;
    
    % Repp 2002
%     The ratio of lag- 0 over lag-1 CCs between ITIs and the IOIs of the pacing signal was computed relative to the lag-1 auto-correlation of the sequence to assess the degree to which individuals predicted upcoming tempo changes (cf. Repp 2002). Prediction indices in these participant samples, which comprised mostly amateur musicians, indicated that about two-thirds of the individuals mainly predicted (ratio[1) ongoing tempo changes. The remaining individuals showed weaker prediction (ratio\1; i.e., they in fact displayed a tendency to track the changes).

%     start_steps     = 2;
%     %start_steps_t   = Time_Steps(start_steps);
%     start_trigs     = find(abs(Time_Trigs-Time_Steps(start_steps))==min(abs(Time_Trigs-Time_Steps(start_steps))));
%     %start_trigs     = find(Time_Trigs<=start_steps_t,1,'last');
%     end_steps       = length(IOI_Steps);
%     end_trigs       = length(IOI_Steps) + start_trigs - start_steps;
%     if doubletrigs == 1
%         IOI_Trigs_temp = IOI_Trigs/2;
%         if mod(length(IOI_Trigs_temp),2)
%             IOI_Trigs_temp=IOI_Trigs_temp(1:end-1);
%         end
%         IOI_Trigs = IOI_Trigs_temp(1:2:end)+IOI_Trigs_temp(2:2:end);
%         l = min([floor(Idx_End_Trigs_IOI/2)-ceil(Idx_Start_Trigs_IOI/2) Idx_End_Steps_IOI-Idx_Start_Steps_IOI]);
%         A = IOI_Trigs(ceil(Idx_Start_Trigs_IOI/2):ceil(Idx_Start_Trigs_IOI/2)+l);
%         M = IOI_Steps(Idx_Start_Steps_IOI:Idx_Start_Steps_IOI+l);
%         A_bl = A - nanmean([A M]);
%         M_bl = M - nanmean([A M]); 
%     else
%         l = min([Idx_End_Trigs_IOI-Idx_Start_Trigs_IOI Idx_End_Steps_IOI-Idx_Start_Steps_IOI]);
%         A = IOI_Trigs(Idx_Start_Trigs_IOI:Idx_Start_Trigs_IOI+l);
%         M = IOI_Steps(Idx_Start_Steps_IOI:Idx_Start_Steps_IOI+l);
%         A_bl = A - nanmean([A M]);
%         M_bl = M - nanmean([A M]);
%     end



%%% NOT GOOD, do it with IOI_nodeviant
%     l = min([Idx_End_Trigs_IOI-Idx_Start_Trigs_IOI Idx_End_Steps_IOI-Idx_Start_Steps_IOI]);
%         A = IOI_Trigs(Idx_Start_Trigs_IOI:Idx_Start_Trigs_IOI+l);
%         M = IOI_Steps(Idx_Start_Steps_IOI:Idx_Start_Steps_IOI+l);
%         A_bl = A - nanmean([A M]);
%         M_bl = M - nanmean([A M]);
%     % IOI_Trigs_nodeviant
%     
%     
%     [r_IOI,lag_IOI] = xcorr(M_bl,A_bl,'coeff');
%     r0_Repp=r_IOI(1,lag_IOI==0);
%     r1_Repp=r_IOI(1,lag_IOI==1);
%     
%     [rac,lagac]=xcorr(A_bl,A_bl,'coeff');
%     ac1_Repp=rac(1,lagac==1);
%     
    
 %% PLOT  
 if fig == 1   

    Line_Trigs              = zeros(size(time));
    Line_Trigs(Idx_Trigs)   = 10;
    Line_Steps              = zeros(size(time));
    Line_Steps(Idx_Steps)   = 10; 
    % plot the data and the generated time-IOI waveforms
    figure; 
    subplot(2,2,1); 
    plot(time,Line_Trigs,'color',[.9 .9 .9]);
    hold on;
    plot(time,Line_Steps,'color',[.8 .8 .8]);
    plot(time,Dots_Trigs,'*b'); plot(time,Dots_Steps,'*r');
    plot(time,Wave_Trigs_keep,'-b');
    plot(time,Wave_Steps_keep,'-r');
    ylim([min([Dots_Trigs Dots_Steps])-.1 max([Dots_Trigs Dots_Steps])+.1])
    
    legend('Targets','Steps','Auditory Targets','Motor Outputs','start analysis','stop analysis');
    start_graph = NaN(1,length(time));
    stop_graph = NaN(1,length(time));
    start_graph(start) = 10; start_graph(start+1) = 0;
    stop_graph(stop) = 10; stop_graph(stop+1) = 0;
    plot(time,start_graph,'g')
    plot(time,stop_graph,'r')
    hold off
    
    subplot(2,2,2); 
    plot(lagtime_wave,corrmax,'color',[.9 .9 .9])
    hold on;
    plot(lagtime_wave,r_wave)
    hold off
    ylim([-1 1]);
    
%     subplot(2,2,3)
%     hold on; plot(A); plot(M); hold off
%     legend('Auditory Targets','Motor Outputs');

    
%     subplot(2,2,4)
%     hold on
%     CorrMax_Repp = lag_IOI(r_IOI==max(r_IOI));
%     corrmax = lag_IOI==CorrMax_Repp;
%     plot(lag_IOI,corrmax,'color',[.9 .9 .9])
%     plot(lag_IOI,r_IOI)
%     ylim([-1 1]);
%     hold off
 end
end

