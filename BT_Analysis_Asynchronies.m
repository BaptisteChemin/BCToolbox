clear; fs                                   = 250; % sampling rate was 250 Hz
% load the .mat file with accelero signals
[matnames, matpath]                         = uigetfile('*.mat', 'Get the stimulation .wav file','MultiSelect','on');
if ischar (matnames); matnames              = {char(matnames)}; end; clc;


for f=1:length(matnames)
    
    matname                                 = char(matnames(f));
    filename                                = [char(matpath) char(matname)];
    dataname                                = matname(1:end-14);
    disp (dataname);
    
    %% get the data
    % use a specific function "GetLatenciesData" to get the data acurately
    % acording to the way they were encoded
    [Time_Steps,Time_Trigs,stimtype]        = XPMCS_DataTransform_GetLatenciesData(filename,'b'); %'b' for both legs, 'l' or 'r' for selected leg
    % select the segment to analyse
    [Time_Steps,Time_Trigs]                 = XPMCS_DataTransform_LatenciesDataCrop(Time_Steps,Time_Trigs,fs); close;
    % time vector
    AllEvents                               = sort([Time_Trigs Time_Steps]);
    time                                    = AllEvents(1)-1:1/fs:AllEvents(end)+1;
    % info about segments to remove due to missing data
    [time_keep] = BT_DataTransform_Time2IEI(Time_Steps,Time_Trigs,time,'removedeviants');
    
    %% MAIN FIGURE
    figure;
    subplot(4,1,1)
    BT_Figure_TimeIOI(Time_Steps,Time_Trigs,time,fs,time_keep)
    % legend(['N aud = ' num2str(length(Times_Trig))],['N steps = ' num2str(length(Times_Step))])
    
    %% ANALYSIS
    %%%%% Probability to observe this time serie; and values of stimulation (ITI)
    if ~strcmp(stimtype,'S')
        [Time_Steps_Surr]                   = BT_Analysis_MonteCarloPhaseRandom(Time_Steps,Time_Trigs,fs,10000,'IOI_type');
        [pA]                                = BT_Analysis_RealDataVsPhaseRandom(Time_Steps,Time_Trigs,Time_Steps_Surr,fs,'MNA',1);  % p MeanAsynch; p stdA; p sumA;
        [pC]                                = BT_Analysis_RealDataVsPhaseRandom(Time_Steps,Time_Trigs,Time_Steps_Surr,fs,'CIRC',0); % p CircStat Rlenght; p tetha;
        legend(['N aud = ' num2str(length(Time_Trigs))],['N steps = ' num2str(length(Time_Steps))],['p obs CS= ' num2str(pC(1,1))], ['p obs SA= ' num2str(pA(1,3))])
        ITI_real                            = BT_Analysis_CircStat_adjust(Time_Trigs,0);
        ITI_mean                            = mean(diff(Time_Trigs));
        ISI_real                            = BT_Analysis_CircStat_adjust(Time_Steps,0);
        ISI_mean                            = mean(diff(Time_Steps));
    else
        if contains(matname,'_T1_')
            ITI_test                        = .712;
        elseif contains(matname,'_T2_')
            ITI_test                        = .619;
        elseif contains(matname,'_T3_') || contains(matname,'_T4_')
            clc; disp(matname(1:2));
            ITI_test                        = input('Spontaneous Gait Rate for participant?   (s) :  ');
        end
        Time_Trigs                          = Time_Trigs(1):ITI_test:Time_Steps(end);
    end
    
    %%%%% Negative Asynchrony
    [MA,stdA,As,Ad]                         = BT_Analysis_NegAsy(Time_Steps,Time_Trigs); % asynchronies, in second
    %%%%% Circular Statistics
    
    [R,rads]                                = BT_Analysis_CircStat(Time_Steps,Time_Trigs,'Mean',0);
    [Rbest,radsbest]                        = BT_Analysis_CircStat(Time_Steps,Time_Trigs,'Best',0);
    
    Rad_Steps_SelfTempo                     = rads(:,1);
    subplot(4,1,2)
    circ_plot(Rad_Steps_SelfTempo,'pretty');
    title('Dispertion around self tempo')
    text(.6,.95,['R = ' num2str(round(R(4,1)*100)/100)])
    text(.75,.8,['ISI = ' num2str(round(R(1,1)*1000)/1000)])
    
    Rad_Steps_TargetTempo   = rads(:,3);
    Rad_Steps_TargetTempo(isnan(Rad_Steps_TargetTempo))=[];
    subplot(4,1,3)
    circ_plot(Rad_Steps_TargetTempo,'pretty');
    title('Dispertion around target')
    text(.6,.95,['R = ' num2str(round(R(4,3)*100)/100)])
    text(.75,.8,['ITI = ' num2str(round(R(1,3)*1000)/1000)])
    
    Index = R(4,3)/R(4,1);
    text(.75,-.95,['Index = ' num2str(Index)])
    
    %%%%% Cross Correlations
    if ~strcmp(stimtype,'S') && ~strcmp(stimtype,'M') && ~strcmp(stimtype,'I')
        % for figure comparison:
        [r_wave,lagtime_wave,rMax_wave_Time,rMax_wave,r0_wave,rac_wave,lagactime_wave] = BT_Analysis_CrossCor(Time_Steps,Time_Trigs,fs,'removedeviants',8,0);
        
        subplot(4,1,4)
        plot(lagactime_wave,rac_wave,'color',[.8 .8 .8])
        hold on
        plot(lagtime_wave,r_wave,'color','b')
        title ('Cross-Correlation function')
        text(rMax_wave_Time,rMax_wave,['\leftarrow r max = ' num2str(rMax_wave) '   lag = ' num2str(rMax_wave_Time)])
        xlim([-10 10])
        
%         Keller      = r0/r1;
%         Repp_tr     = (r1-ac1)/(1-ac1);
%         Repp_pr     = (r0-ac1)/(1-ac1);
        %%%%% ACCorr = (r_wave+1)./(rac_wave+1); %%%%%%

        %[p] = BT_Analysis_RealDataVsPhaseRandom(Time_Steps,Time_Trigs,Time_Steps_Surr,fs,'CROSSCO',1);
    end
    %%%
    
    %%%%% frequency
    central_IOI = BT_Analysis_CircStat_adjust(Time_Steps,0);
    Frequency_ratio = central_IOI/(ITI*2);
    
    
    %%%
    
    
    %% Save the RESULTS
    % Figure
    ID = [char(fichiers(file,1)) '_' char(fichiers(file,2))];
    savefig(['XPPARKGAIT_RESULTS_Fig_' char(ID) '.fig']);
    % Values
    
    RESULTS.Asynchronies.MNA    = MA;
    RESULTS.Asynchronies.stdNA  = stdA;
    RESULTS.Asynchronies.NAs    = As;
    % RESULTS.Asynchronies.pvalue = p(1,3);
    
    RESULTS.CircStats.R         = R;
    RESULTS.CircStats.rads      = rads;
    RESULTS.CircStats.Index     = Index;
    
    if ~strcmp(title_maker(5:7),'Spo')  && ~strcmp(title_maker(5:7),'Hap') && ~strcmp(title_maker(5:7),'Iso')
        RESULTS.CrossCorr.lagtime_wave  = lagtime_wave;
        RESULTS.CrossCorr.r_wave        = r_wave;
        RESULTS.CrossCorr.acMax_wave    = acMax_wave;
        RESULTS.CrossCorr.r0            = r0;
        RESULTS.CrossCorr.r1            = r1;
        RESULTS.CrossCorr.ac1           = ac1;
    end
    
    RESULTS.FrequencyRatio.ISI          = central_IOI;
    RESULTS.FrequencyRatio.ITI          = ITI;
    RESULTS.FrequencyRatio.ratio        = Frequency_ratio;
    
    if ~strcmp(title_maker(5:7),'Spo')
        RESULTS.surrogates              = Time_Steps_Surr;
    end
    
    eval(sprintf('%s',['save ([''XPPARKGAIT_RESULTS_Data' char(ID) '.mat''], ''RESULTS'',''ID'');']));
    close
end
