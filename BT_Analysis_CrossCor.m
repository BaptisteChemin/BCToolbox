function [r_wave,lagtime_wave,rMax_wave_Time,rMax_wave,r0_wave,acMax_wave,r0_Repp,r1_Repp,ac1_Repp] = BT_Analysis_CrossCor(Time_Steps,Time_Trigs,fs)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    %% data
    time                    = (Time_Trigs(1)-1):1/fs:(Time_Trigs(end)+1);
    
    IOI_Trigs               = diff(Time_Trigs);
    IOI_Steps               = diff(Time_Steps);

    [idx_steps,idx_trigs]   = BT_DataTransform_TimeIdx(Time_Steps,Time_Trigs,time);
    
    Line_Trigs              = zeros(size(time));
    Line_Trigs(idx_trigs)   = 10;
    Line_Steps              = zeros(size(time));
    Line_Steps(idx_steps)   = 10;
    
    %% limits of analysis
    Idx_Start_Steps_time    = 2; % Index of the second actual step
    %Time_Start_Steps        = Time_Steps(Idx_Start_Steps_time); % latency of the second actual step
    Idx_Start_Trigs_time    = find(abs(Time_Trigs-Time_Steps(Idx_Start_Steps_time))==min(abs(Time_Trigs-Time_Steps(Idx_Start_Steps_time)))); % Index of the closest event to second step
    %Time_Start_Trigs        = Time_Trigs(Idx_Start_Steps_time); % latency of the second closest event to second step
    
    Idx_End_Trigs_time      = length(Time_Trigs); % Index of the last stimulus
    %Time_End_Trigs          = Time_Trigs(Idx_End_Trigs_time); % latency of the last stimulus
    Idx_End_Steps_time      = find(abs(Time_Steps-Time_Trigs(Idx_End_Trigs_time))==min(abs(Time_Steps-Time_Trigs(Idx_End_Trigs_time)))); % Index of the closest steps to last stimulus
    %Time_End_Steps          = Time_Steps(Idx_End_Steps_time); % latency of the closest steps to last stimulus
    
    Idx_Start_Steps_IOI     = Idx_Start_Steps_time - 1; if Idx_Start_Steps_IOI==0; Idx_Start_Steps_IOI=1; end
    Idx_Start_Trigs_IOI     = Idx_Start_Trigs_time - 1; if Idx_Start_Trigs_IOI==0; Idx_Start_Trigs_IOI=1; end
    Idx_End_Steps_IOI       = Idx_End_Steps_time - 1;
    Idx_End_Trigs_IOI       = Idx_End_Trigs_time - 1;
    
    % check the limits? 
    
    %% cross-correlations on time-IOI waveforms
    % generate the waveforms
    Dots_Trigs = NaN(1,length(time));
    for f=1:length(IOI_Trigs)
        Dots_Trigs(idx_trigs(f+1))=IOI_Trigs(f); 
    end
    Wave_Trigs = fillmissing(Dots_Trigs,'linear');
    
    Dots_Steps = NaN(1,length(time));
    for f=1:length(IOI_Steps)
        Dots_Steps(idx_steps(f+1))=IOI_Steps(f);
    end
    % deal with "missing events" by supressing the large IOI
    outliers = isoutlier(Dots_Steps,'gesd');
    Dots_Steps_corrected = Dots_Steps;
    Dots_Steps_corrected(outliers==1) = NaN; 
    
    Wave_Steps = fillmissing(Dots_Steps_corrected,'linear');
    % crop the waveforms so there is no extrapolation due to fillmissing function
%     Wave_Trigs(1:idx_trigs(2)-.3*fs) = NaN;
%     Wave_Trigs(idx_trigs(end)+.3*fs:end) = NaN;
%    
%     Wave_Steps(1:idx_steps(2)-.3*fs)=NaN; 
%     Wave_Steps(idx_steps(end)+.3*fs:end)=NaN;
    
    % plot the data and the generated time-IOI waveforms
    figure; 
    subplot(2,2,1); 
    plot(time,Line_Trigs,'color',[.9 .9 .9]);
    hold on;
    plot(time,Line_Steps,'color',[.8 .8 .8]);
    plot(time,Dots_Trigs,'*b'); plot(time,Dots_Steps,'*r');
    plot(time,Wave_Trigs,'-b');
    plot(time,Wave_Steps,'-r');
    ylim([min([Dots_Trigs Dots_Steps])-.1 max([Dots_Trigs Dots_Steps])+.1])
    
    legend('Targets','Steps','Auditory Targets','Motor Outputs','start analysis','stop analysis');
    
    % correlations of time-IOI waveforms, and maximum correlation lag in
    % time.
    start       = idx_steps(Idx_Start_Steps_time)+1;
    stop        = idx_steps(Idx_End_Steps_time)-1;
    
    start_graph = NaN(1,length(time));
    stop_graph = NaN(1,length(time));
    start_graph(start) = 10; start_graph(start+1) = 0;
    stop_graph(stop) = 10; stop_graph(stop+1) = 0;
    plot(time,start_graph,'g')
    plot(time,stop_graph,'r')
    hold off
        
    Wave_Trigs_bl = Wave_Trigs - nanmean([Wave_Trigs Wave_Steps]);
    Wave_Steps_bl = Wave_Steps - nanmean([Wave_Trigs Wave_Steps]);
    
    [r_wave,lag_wave]=xcorr(Wave_Steps_bl(start:stop),Wave_Trigs_bl(start:stop),'coeff');
    lagtime_wave = lag_wave/fs;
    
    window = [find(lagtime_wave==-20) find(lagtime_wave==+20)];
    rMax_wave_Time = lagtime_wave(find(r_wave(window(1):window(2))==max(r_wave(window(1):window(2))))+ window(1));
    rMax_wave = max(r_wave(window(1):window(2)));
    corrmax = lagtime_wave==rMax_wave_Time;
    
    subplot(2,2,2); 
    plot(lagtime_wave,corrmax,'color',[.9 .9 .9])
    hold on;
    plot(lagtime_wave,r_wave)
    hold off
    ylim([-1 1]);
    
    r0_idx = lagtime_wave==0;
    r0_wave = r_wave(r0_idx);
    
    [rac_wave,lagac_wave]=xcorr(Wave_Trigs_bl(start:stop),Wave_Trigs_bl(start:stop),'coeff');
    lagactime_wave = lagac_wave/fs;
    
    acMax_wave=rac_wave(1,lagactime_wave==rMax_wave_Time);
    
    % Repp 2002
%     The ratio of lag- 0 over lag-1 CCs between ITIs and the IOIs of the pacing signal was computed relative to the lag-1 auto-correlation of the sequence to assess the degree to which individuals predicted upcoming tempo changes (cf. Repp 2002). Prediction indices in these participant samples, which comprised mostly amateur musicians, indicated that about two-thirds of the individuals mainly predicted (ratio[1) ongoing tempo changes. The remaining individuals showed weaker prediction (ratio\1; i.e., they in fact displayed a tendency to track the changes).

%     start_steps     = 2;
%     %start_steps_t   = Time_Steps(start_steps);
%     start_trigs     = find(abs(Time_Trigs-Time_Steps(start_steps))==min(abs(Time_Trigs-Time_Steps(start_steps))));
%     %start_trigs     = find(Time_Trigs<=start_steps_t,1,'last');
%     end_steps       = length(IOI_Steps);
%     end_trigs       = length(IOI_Steps) + start_trigs - start_steps;
    l = min([Idx_End_Trigs_IOI-Idx_Start_Trigs_IOI Idx_End_Steps_IOI-Idx_Start_Steps_IOI]);
    A = IOI_Trigs(Idx_Start_Trigs_IOI:Idx_Start_Trigs_IOI+l);
    M = IOI_Steps(Idx_Start_Steps_IOI:Idx_Start_Steps_IOI+l);
    
    subplot(2,2,3)
    hold on; plot(A); plot(M); hold off
    legend('Auditory Targets','Motor Outputs');
    A_bl = A - nanmean([A M]);
    M_bl = M - nanmean([A M]);
    
    [r_IOI,lag_IOI] = xcorr(M_bl,A_bl,'coeff');
    r0_Repp=r_IOI(1,lag_IOI==0);
    r1_Repp=r_IOI(1,lag_IOI==1);
    
    [rac,lagac]=xcorr(A_bl,A_bl,'coeff');
    ac1_Repp=rac(1,lagac==1);
    
    subplot(2,2,4)
    hold on
    CorrMax_Repp = lag_IOI(r_IOI==max(r_IOI));
    corrmax = lag_IOI==CorrMax_Repp;
    plot(lag_IOI,corrmax,'color',[.9 .9 .9])
    plot(lag_IOI,r_IOI)
    ylim([-1 1]);
    hold off
end

