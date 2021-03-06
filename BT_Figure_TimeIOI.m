function BT_Figure_TimeIOI(Time_Steps,Time_Trigs,time,fs,time_keep)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    

    f_avail = exist('yyaxis','file');
    %% make sure the data are suited for the graph
    Time_Steps(Time_Steps<=time(2)) = [];
    if isempty(Time_Trigs)
        Time_Trigs = [time(2) time(end-1)];
    end
    
    %% plot each event as a line at its time latency
    [idx_steps,idx_trigs] = BT_DataTransform_Time2Idx(Time_Steps,Time_Trigs,time);
    Line_Trigs      = zeros(size(time));
    Line_Trigs(idx_trigs) = 10;
    Line_Steps      = zeros(size(time));
    Line_Steps(idx_steps) = 10;
    
    if f_avail~=0; yyaxis left; end
    plot(time,Line_Trigs,'color',[.95 .525 .298],'LineStyle','-');
    hold on;
    plot(time,Line_Steps,'color',[.4 .677 .971],'LineStyle','-');
    
    
    %% plot the relative phase
    % compute the relative phase
    if f_avail~=0
        Deg_Rel_TrigSteps           = [];
        for s=1:length(Time_Steps)
            Trig_Before_Step        = Time_Trigs(find(Time_Trigs<Time_Steps(s),1,'last'));
            Trig_After_Step         = Time_Trigs(find(Time_Trigs>=Time_Steps(s),1,'first'));
            if ~isempty(Trig_Before_Step) && ~isempty(Trig_After_Step)
                Instant_Period      = Trig_After_Step-Trig_Before_Step;
                Instant_Phase       = circ_rad2ang(2*pi*((Time_Steps(s)-Trig_Before_Step)/Instant_Period));
                if 180 <= Instant_Phase && Instant_Phase < 360 % keep instant phase within a -180 to 180 � range
                    Instant_Phase = Instant_Phase-360;
                end
                Deg_Rel_TrigSteps   = [Deg_Rel_TrigSteps Instant_Phase]; %#ok<*AGROW>
                %Rad_Rel_TrigSteps   = [Deg_Rel_TrigSteps 2*pi*((Time_Steps(s)-Trig_Before_Step)/Instant_Period)]; %#ok<*AGROW>
            else
                Deg_Rel_TrigSteps   = [Deg_Rel_TrigSteps NaN];
            end
        end
        Dots_Phase = NaN(1,length(time));
        for f=1:length(Deg_Rel_TrigSteps)
            Dots_Phase(idx_steps(f))=Deg_Rel_TrigSteps(f);
        end
        
        yyaxis right
        plot(time,Dots_Phase,'s','color',[.5 .5 .5]); hold on;
        plot(time,fillmissing(Dots_Phase,'linear'),'color',[.5 .5 .5]);
        ylim([-180 1000]);
        yticks([-180 0 180]);
    end
    
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

    if f_avail~=0; yyaxis left; end
    plot(time,Dots_Trigs,'diamond','color',[.85 .325 .098],'LineWidth',2); plot(time,Dots_Steps,'o','color',[0 .447 .741],'LineWidth',2);
    
    %% time-IOI waveforms
    if f_avail~=0
        if ~exist('time_keep','var')
            time_keep = time;
        end

        Wave_Steps = fillmissing(Dots_Steps,'linear');
        Wave_Trigs = fillmissing(Dots_Trigs,'linear');
        % crop the waveforms so there is no extrapolation due to fillmissing function
        Wave_Trigs(1:idx_trigs(2)-.3*fs) = NaN;
        Wave_Trigs(idx_trigs(end)+.3*fs:end) = NaN;
        
        Wave_Steps(1:idx_steps(2)-.3*fs)=NaN;
        Wave_Steps(idx_steps(end)+.3*fs:end)=NaN;
        
        % remove data that are not kept in time_keep vector
        Wave_Steps(isnan(time_keep)) = NaN;
        Wave_Trigs(isnan(time_keep)) = NaN;
        
        % plot the time-IOI waveforms
        yyaxis left
        plot(time,Wave_Trigs,'color',[.85 .325 .098],'LineWidth',2,'LineStyle','-');
        plot(time,Wave_Steps,'color',[0 .447 .741],'LineWidth',2,'LineStyle','-');
    end
%     if length(Time_Trigs)>2
%         ylim([min([Dots_Trigs Dots_Steps])-.1 max([Dots_Trigs Dots_Steps])+.1])
%     else
%         ylim([min(Dots_Steps)-.1 max(Dots_Steps)+.1])
%     end
    ylim([0 2]);
    xlim([time(1) time(end)])
    
end

