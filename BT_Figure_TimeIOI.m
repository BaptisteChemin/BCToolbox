function BT_Figure_TimeIOI(Time_Steps,Time_Trigs,time,fs)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    
    %time            = (Time_Trigs(1)-1):1/fs:(Time_Trigs(end)+1);
    %% plot each event as a line at its time latency
    Time_Steps(Time_Steps<=time(2)) = [];
    
    [idx_steps,idx_trigs] = BT_DataTransform_TimeIdx(Time_Steps,Time_Trigs,time);
    Line_Trigs      = zeros(size(time));
    Line_Trigs(idx_trigs) = 10;
    Line_Steps      = zeros(size(time));
    Line_Steps(idx_steps) = 10;
    
    plot(time,Line_Trigs,'color',[.95 .525 .298]);
    hold on;
    plot(time,Line_Steps,'color',[.4 .677 .971]);
    
    %% plot each event as a dot at its time-latency with a value = the preceding InterEventInterval
    IOI_Trigs = diff(Time_Trigs);
    IOI_Steps = diff(Time_Steps);
    
    Dots_Trigs = NaN(1,length(time));
    for f=1:length(IOI_Trigs)
        Dots_Trigs(idx_trigs(f+1))=IOI_Trigs(f);  
    end

    Dots_Steps = NaN(1,length(time));
    for f=1:length(IOI_Steps)
        Dots_Steps(idx_steps(f+1))=IOI_Steps(f);
    end

    plot(time,Dots_Trigs,'diamond','color',[.85 .325 .098],'LineWidth',2); plot(time,Dots_Steps,'o','color',[0 .447 .741],'LineWidth',2);
    
    %% time-IOI waveforms
    % deal with "missing events" by supressing the large IOI
    outliers = isoutlier(Dots_Steps,'gesd');
    Dots_Steps_corrected = Dots_Steps;
    Dots_Steps_corrected(outliers==1) = NaN;
    Wave_Steps = fillmissing(Dots_Steps_corrected,'linear');
    Wave_Trigs = fillmissing(Dots_Trigs,'linear');
    % crop the waveforms so there is no extrapolation due to fillmissing function
    
    Wave_Trigs(1:idx_trigs(2)-.3*fs) = NaN;
    Wave_Trigs(idx_trigs(end)+.3*fs:end) = NaN;
   
    Wave_Steps(1:idx_steps(2)-.3*fs)=NaN; 
    Wave_Steps(idx_steps(end)+.3*fs:end)=NaN;
    
    % plot the time-IOI waveforms
    plot(time,Wave_Trigs,'color',[.85 .325 .098],'LineWidth',2);
    plot(time,Wave_Steps,'color',[0 .447 .741],'LineWidth',2);
    
    ylim([min([Dots_Trigs Dots_Steps])-.1 max([Dots_Trigs Dots_Steps])+.1])
end

