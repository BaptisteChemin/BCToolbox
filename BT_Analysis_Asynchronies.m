clear; 
fs                                      = 250; % sampling rate was 250 Hz
% load the .mat file with accelero signals
[matnames, matpath]                     = uigetfile('*.mat', 'Get the stimulation .wav file','MultiSelect','on');
if ischar (matnames); matnames          = {char(matnames)}; end; clc;

%% First, take the Spontaneous Gait time serie and keep it as reference
S = find(contains(matnames,'_S_')==1);
matname                                 = char(matnames(S));
filename                                = [char(matpath) char(matname)];
dataname                                = matname(1:end-14);
[Time_Steps_S,Time_Trigs_S,stimtype]    = XPMCS_DataTransform_GetLatenciesData(filename,'b'); %'b' for both legs, 'l' or 'r' for selected leg
% select the segment to analyse
[Time_Steps_S,Time_Trigs_S]             = XPMCS_DataTransform_LatenciesDataCrop(Time_Steps_S,Time_Trigs_S,fs); close;
% time vector
AllEvents_S                             = sort([Time_Trigs_S Time_Steps_S]);
time_S                                  = AllEvents_S(1)-1:1/fs:AllEvents_S(end)+1;
% info about segments to remove due to missing data
[time_keep_S]                           = BT_DataTransform_Time2IEI(Time_Steps_S,Time_Trigs_S,time_S,'removedeviants');
figure; BT_Figure_TimeIOI               (Time_Steps_S,Time_Trigs_S,time_S,fs,time_keep_S)   
% generate surrogate time series with the Spontaneous time serie
[Time_Steps_S_Surr]                     = BT_Analysis_MonteCarloPhaseRandom(Time_Steps_S,Time_Trigs,fs,10000,'wave_type');

%% Loop the analysis (and figure generation) for each trial
for trial = 1:length(matnames)
    matname                           	= char(matnames(trial));
    filename                            = [char(matpath) char(matname)];
    dataname                            = matname(1:end-14);
    disp (dataname);
    
    %% get the data
    % Use a specific function "GetLatenciesData" to properly retrieve the data, acording to the way they were encoded
    [Time_Steps,Time_Trigs,stimtype]    = XPMCS_DataTransform_GetLatenciesData(filename,'b'); %'b' for both legs, 'l' or 'r' for selected leg
    % Select the segment to analyse
    [Time_Steps,Time_Trigs]             = XPMCS_DataTransform_LatenciesDataCrop(Time_Steps,Time_Trigs,fs); close;
    % Create a time vector
    AllEvents                           = sort([Time_Trigs Time_Steps]);
    time                                = AllEvents(1)-1:1/fs:AllEvents(end)+1;
    % Get info about segments to remove due to missing values. 
    % Note that all the following analysis are able to deal with sporadic missing values (not tested above 3 consecutive missing events)
    [time_keep]                         = BT_DataTransform_Time2IEI(Time_Steps,Time_Trigs,time,'removedeviants');
    
    %% MAIN FIGURE
    figure;
    %%%%% Time-Frequency Representation
    subplot(4,1,1)
    BT_Figure_TimeIOI                   (Time_Steps,Time_Trigs,time,fs,time_keep);    
    legend                              (['N aud = ' num2str(length(Time_Trigs))],['N steps = ' num2str(length(Time_Steps))])
    %%%%% Circular Statistics
    [R,rads]                            = BT_Analysis_CircStat(Time_Steps,Time_Trigs,'Median',0); %     [R,rads] = BT_Analysis_CircStat(Time_Steps,Time_Trigs,'Mean',1);
    
    Rad_Steps_SelfTempo                 = rads(:,1);
    subplot(4,1,2)
    circ_plot                           (Rad_Steps_SelfTempo,'pretty');
    title                               ('Dispertion around self tempo')
    text                                (.6,.95,['R = ' num2str(round(R(4,1)*100)/100)])
    text                                (.75,.8,['ISI = ' num2str(round(R(1,1)*1000)/1000)])
    
    Rad_Steps_TargetTempo               = rads(:,3);
    Rad_Steps_TargetTempo               (isnan(Rad_Steps_TargetTempo))=[];
    subplot(4,1,3)
    circ_plot                           (Rad_Steps_TargetTempo,'pretty');
    title                               ('Dispertion around target')
    text                                (.6,.95,['R = ' num2str(round(R(4,3)*100)/100)])
    text                                (.75,.8,['ITI = ' num2str(round(R(1,3)*1000)/1000)])
    
    Index                               = R(4,3)/R(4,1);
    text                                (.75,-.95,['Index = ' num2str(Index)])
    
    %%%%% Cross Correlations
    if ~strcmp(stimtype,'S') && ~strcmp(stimtype,'I') % && ~strcmp(stimtype,'M')
        [r_wave,lagtime_wave,rMax_wave_Time,rMax_wave,r0_wave,rac_wave,lagactime_wave] = BT_Analysis_CrossCor(Time_Steps,Time_Trigs,fs,'removedeviants',5,0);
        subplot(4,1,4); hold on
        plot                            (lagactime_wave,rac_wave,'color',[.8 .8 .8])
        plot                            (lagtime_wave,r_wave,'color','b')
        title                           ('Cross-Correlation function')
        text                            (rMax_wave_Time,rMax_wave,['\leftarrow r max = ' num2str(rMax_wave) '   lag = ' num2str(rMax_wave_Time)])
        xlim                            ([-20 20])          
%     [r,lagtime,~,~,~,rac,lagactime]   = BT_Analysis_CrossCor(Time_Steps,Time_Trigs,fs,time_keep,5,0);
%     figure; plot(lagtime,r,'DisplayName','r_wave_obs');hold on;plot(lagactime,rac,'DisplayName','rac_wave');hold off;
%     % When tracking r1 is maximal, r0 will approach ac1 because ITIs echo
%     % the IOIs at a lag of 1.
%     % Therefore, a correction is applied to both r0 and r1,
%     % in order to take into account the fact that different timing patterns
%     % have different ac1 values.
%      r0corr = (r0-ac1)/(1-ac1)
%      r1corr = (r1-ac1)/(1-ac1) 
    end
    %%%%% Negative Asynchrony
    [MA,stdA,As,Ad]                     = BT_Analysis_NegAsy(Time_Steps,Time_Trigs); % asynchronies, in second
    
    %% ANALYSIS
    %%%%% Test for Phase Alignment
    if ~strcmp(stimtype,'S')
        [Time_Steps_Surr]               = BT_Analysis_MonteCarloPhaseRandom(Time_Steps,Time_Trigs,fs,10000,'wave_type');
        [pA]                            = BT_Analysis_RealDataVsPhaseRandom(Time_Steps,Time_Trigs,Time_Steps_Surr,fs,'MNA',0);  % p MeanAsynch; p stdA; p sumA;
        [pC]                            = BT_Analysis_RealDataVsPhaseRandom(Time_Steps,Time_Trigs,Time_Steps_Surr,fs,'CIRC',0); % p CircStat Rlenght; p tetha;
        [pSI]                           = BT_Analysis_RealDataVsPhaseRandom(Time_Steps,Time_Trigs,Time_Steps_Surr,fs,'SI',0); % p Synchr Index
        [pXC]                           = BT_Analysis_RealDataVsPhaseRandom(Time_Steps,Time_Trigs,Time_Steps_Surr,fs,'XCORR',0); % p(r), p(lag)
    end
    
    %%%%% Test for Frequency Alignement      
    % Test the distribution of instantaneous periods of S gait VS the distribution of instantaneous periods of paced gait.
    [~,~,IOI_Spont]                     = BT_DataTransform_Time2IEI(Time_Steps_S,Time_Trigs,time,'removedeviants');
    [~,IOI_Trigs,IOI_Paced]             = BT_DataTransform_Time2IEI(Time_Steps  ,Time_Trigs,time,'removedeviants');
    [pW,hW,statsW]                      = ranksum(IOI_Spont,IOI_Paced);
    IOI_Spont_median                    = median(IOI_Spont);
    IOI_Paced_median                    = median(IOI_Paced);
    IOI_Trigs_median                    = median(IOI_Trigs);
    
%     % Pre-Visualisation with frequency ratios
%     if ~strcmp(stimtype,'S')
%         ITI_best                        = BT_Analysis_CircStat_adjust(Time_Trigs,0);
%         ITI_mean                        = mean(diff(Time_Trigs));
%         ITI_median                      = median(diff(Time_Trigs));
%     end
%     ISI_best                            = BT_Analysis_CircStat_adjust(Time_Steps,0);
%     ISI_mean                            = mean(diff(Time_Steps));
%     ISI_median                          = median(diff(Time_Steps));
%     
%     ISIsp_best                          = BT_Analysis_CircStat_adjust(Time_Steps_S,0);
%     ISIsp_mean                          = mean(diff(Time_Steps_S));
%     ISIsp_median                        = median(diff(Time_Steps_S));
%     
%     Frequency_Ratio_target              = ISI_median/ITI_median; % This ratio is close to 1 if the gait frequency is well alignend to the beat frequency
%     Frequency_Ratio_spontaneous         = ISIsp_median/ISI_median; % This ratio is close to 1 if the spontaneous gait frequency is similar to the cued gait frequency
%     Adjustment_Factor                   = 100*(Frequency_Ratio_target-Frequency_Ratio_spontaneous)/Frequency_Ratio_spontaneous;
%     
%     % Test of hypothesis "Spontaneous Gait is not synchronized with the
%     % acoustic Stimulus". If pC_S < .05, then the S gait is sync with acoustic stimulus (not good).
%     [R_S,rads_S]                        = BT_Analysis_CircStat(Time_Steps_S,Time_Trigs,'Median',1);
%     [pC_S]                              = BT_Analysis_RealDataVsPhaseRandom(Time_Steps_S,Time_Trigs,Time_Steps_S_Surr,fs,'CIRC',1);
%     
    % Test the distribution of Relative Phase between S gait and acoustic
    % stim VS the distribution of Relative Phase between paced gait and
    % acoustic stim.
    
    
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
