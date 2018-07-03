clear; 
fs                                      = 250; % sampling rate was 250 Hz
% load the .mat file with accelero signals
[matnames, matpath]                     = uigetfile('*.mat', 'Get the stimulation .wav file','MultiSelect','on');
if ischar (matnames); matnames          = {char(matnames)}; end; clc;

%% First, take the Spontaneous Gait time serie and keep it as reference
S = find(contains(matnames,'_S_')==1,1,'first');
matname                                 = char(matnames(S));
filename                                = [char(matpath) char(matname)];
dataname                                = matname(1:end-14);
[Time_Steps_S,Time_Trigs_S,stimtype]    = XPMCS_DataTransform_GetLatenciesData(filename,'b'); %'b' for both legs, 'l' or 'r' for selected leg
% select the segment to analyse
%[Time_Steps_S,Time_Trigs_S]             = XPMCS_DataTransform_LatenciesDataCrop(Time_Steps_S,Time_Trigs_S,fs); close;
[Time_Steps_S,Time_Trigs_S]             = BT_DataTransform_LatenciesDataCropAuto(Time_Steps_S,Time_Trigs_S,fs,2,79);
% time vector
AllEvents_S                             = sort([Time_Trigs_S Time_Steps_S]);
time_S                                  = AllEvents_S(1)-1:1/fs:AllEvents_S(end)+1;
% info about segments to remove due to missing data
[time_keep_S]                           = BT_DataTransform_Time2IEI(Time_Steps_S,Time_Trigs_S,time_S,'removedeviants');
% figure; BT_Figure_TimeIOI               (Time_Steps_S,Time_Trigs_S,time_S,fs,time_keep_S)   
% generate surrogate time series with the Spontaneous time serie
[Time_Steps_S_Surr]                     = BT_Analysis_MonteCarloPhaseRandom(Time_Steps_S,Time_Trigs_S,fs,10000,'wave_type');

%% Second, identify the pacer tempo
I = find(contains(matnames,'_I_')==1);
if ~isempty(I)
    matname                                 = char(matnames(I));
    filename                                = [char(matpath) char(matname)];
    dataname                                = matname(1:end-14);
    [~,Time_Trigs_M]                        = XPMCS_DataTransform_GetLatenciesData(filename,'b');
    stimperiod                              = mean(diff(Time_Trigs_M));
else
    stimperiod                              = str2double(input('What is the Pacer Tempo? (ms):   ','s'))/1000;
end

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
    %[Time_Steps,Time_Trigs]             = XPMCS_DataTransform_LatenciesDataCrop(Time_Steps,Time_Trigs,fs); close;
    [Time_Steps,Time_Trigs]             = BT_DataTransform_LatenciesDataCropAuto(Time_Steps,Time_Trigs,fs,2,79);
    % Create a time vector
    AllEvents                           = sort([Time_Trigs Time_Steps]);
    time                                = AllEvents(1)-1:1/fs:AllEvents(end)+1;
    % Get info about segments to remove due to missing values. 
    % Note that all the following analysis are able to deal with sporadic missing values (not tested above 3 consecutive missing events)
    [time_keep]                         = BT_DataTransform_Time2IEI(Time_Steps,Time_Trigs,time,'removedeviants');
    if ~isempty(find(isnan(time_keep), 1))
         figure;
         BT_Figure_TimeIOI                   (Time_Steps,Time_Trigs,time,fs,time_keep);   
         k = input('discard the events?  (y - n):    ','s'); 
         if strcmp(k,'n')
             time_keep(isnan(time_keep)) = 1;
         end
         close
    end
    %% MAIN FIGURE
    figure;
    %%%%% Time-Frequency Representation
    subplot(4,1,1)
    BT_Figure_TimeIOI                   (Time_Steps,Time_Trigs,time,fs,time_keep);    
    legend                              (['N aud = ' num2str(length(Time_Trigs))],['N steps = ' num2str(length(Time_Steps))])
    %%%%% Circular Statistics
    if strcmp(stimtype,'S')
        Time_Trigs = 0:stimperiod:Time_Steps(end)+stimperiod;
    end
    
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
        [r_wave,lagtime_wave,rMax_wave_Time,rMax_wave,~,rac_wave,lagactime_wave] = BT_Analysis_CrossCor(Time_Steps,Time_Trigs,fs,'removedeviants',5,0);
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
        %[pSI]                           = BT_Analysis_RealDataVsPhaseRandom(Time_Steps,Time_Trigs,Time_Steps_Surr,fs,'SI',0); % p Synchr Index
        if ~strcmp(stimtype,'I')
            [pXC]                           = BT_Analysis_RealDataVsPhaseRandom(Time_Steps,Time_Trigs,Time_Steps_Surr,fs,'XCORR',0); % p(r), p(lag)
        end
    end
    
    %%%%% Test for Frequency Alignement      
    % Test the distribution of instantaneous periods of S gait VS the distribution of instantaneous periods of paced gait.
    [~,~,IOI_Spont]                     = BT_DataTransform_Time2IEI(Time_Steps_S,Time_Trigs,time,'removedeviants');
    %IOI_Spont = diff(Time_Steps_S); if the line above generates a bug
    [~,IOI_Trigs,IOI_Paced]             = BT_DataTransform_Time2IEI(Time_Steps  ,Time_Trigs,time,'removedeviants');
    [pW,~,statsW]                       = ranksum(IOI_Spont,IOI_Paced);
    IOI_Spont_median                    = median(IOI_Spont);
    IOI_Paced_median                    = median(IOI_Paced);
    IOI_Trigs_median                    = median(IOI_Trigs);
         
    % Test of hypothesis "Spontaneous Gait is not synchronized with the
    % acoustic Stimulus". If pC_S < .05, then the S gait is sync with acoustic stimulus 
    % (POSSIBILITY OF FALSE POSITIVE RESULT FOR SYNCHRONIZATION).
    if ~strcmp(stimtype,'S')
        [R_S,rads_S]                        = BT_Analysis_CircStat(Time_Steps_S,Time_Trigs,'Median',0);
        [pC_S]                              = BT_Analysis_RealDataVsPhaseRandom(Time_Steps_S,Time_Trigs,Time_Steps_S_Surr,fs,'CIRC',0);
    end
    % Test the distribution of Relative Phase between I, AC, and NC gait
    % (maybe do it in a separate script, because I'll have all the necessay
    % data)
    
    %% Save the RESULTS
    % Figure
    ID = [char(dataname) '_results'];
    savefig([char(ID) '.fig']);
    % Values
    if ~strcmp(stimtype,'S')
        RESULTS.Asynchronies.MA     = MA;
        RESULTS.Asynchronies.stdA   = stdA;
        RESULTS.Asynchronies.As     = As;
        RESULTS.Asynchronies.Ad     = Ad;
        RESULTS.Asynchronies.pMA    = pA(1,1); %p MeanAsynch; p stdA; p sumA;
        RESULTS.Asynchronies.pstdA  = pA(1,2);
        RESULTS.Asynchronies.psumA  = pA(1,3);
    end
    
    RESULTS.CircStats.Rself     = R(:,1);
    RESULTS.CircStats.Rbeat     = R(:,3);
    if ~strcmp(stimtype,'S')
        RESULTS.CircStats.Rspont    = R_S(:,3);
        RESULTS.CircStats.RadsSelf  = rads(:,1);
        RESULTS.CircStats.RadsBeat  = rads(:,3);
        RESULTS.CircStats.RadsSpont = rads_S(:,3);
        RESULTS.CircStats.Index     = Index;

        RESULTS.CircStats.pRbeat    = pC(1,1);
        RESULTS.CircStats.pTbeat    = pC(1,2);
        RESULTS.CircStats.pRspont   = pC_S(1,1);
        RESULTS.CircStats.pTspont   = pC_S(1,2);
    end
    
    if ~strcmp(stimtype,'S') && ~strcmp(stimtype,'I')
        RESULTS.CrossCorr.lag                   = rMax_wave_Time;
        RESULTS.CrossCorr.rmax                  = rMax_wave;
        RESULTS.CrossCorr.xcorr_functionwav     = r_wave;
        RESULTS.CrossCorr.xcorr_functionlag     = lagtime_wave;
        RESULTS.CrossCorr.ac_functionwav        = rac_wave;
        RESULTS.CrossCorr.ac_functionlag        = lagactime_wave;
        RESULTS.CrossCorr.pxcorrmax             = pXC(1,1);
        RESULTS.CrossCorr.plagmax               = NaN;
    end
    
    RESULTS.FrequencyTest.IOI_Spont             = IOI_Spont;
    RESULTS.FrequencyTest.IOI_Paced             = IOI_Paced;
    RESULTS.FrequencyTest.IOI_Trigs             = IOI_Trigs;
    RESULTS.FrequencyTest.IOI_Spont             = IOI_Spont;
    RESULTS.FrequencyTest.pWilcoxon             = pW;
    RESULTS.FrequencyTest.statsWilcoxon         = statsW; 
    
    RESULTS.latencies.Time_Steps                = Time_Steps;
    RESULTS.latencies.Time_Trigs                = Time_Trigs;
    RESULTS.latencies.Time_Steps_S              = Time_Steps_S;
    if ~strcmp(stimtype,'S')
        RESULTS.latencies.surrogates_syncstep   = Time_Steps_Surr;
        RESULTS.latencies.surrogates_spontstep  = Time_Steps_S_Surr;
    end
    
    eval(sprintf('%s',['save ([''' char(ID) '.mat''], ''RESULTS'',''dataname'');']));
    close all
end
