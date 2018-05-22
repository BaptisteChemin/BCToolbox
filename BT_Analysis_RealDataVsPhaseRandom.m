function [p] = BT_Analysis_RealDataVsPhaseRandom(Time_Steps_Obs,Time_Trigs_Obs,Time_Steps_Surr,fs,type,fig)
%XPGAIT_RealDataVsPhaseRandom test a time serie vs surrogate Data



if strcmp(char(type),'MNA')
    % Negative Asynch
    [MA_obs,stdA_obs,As] = BT_Analysis_NegAsy(Time_Steps_Obs,Time_Trigs_Obs);
    sumNA_obs = nansum(abs(As));
    MA_sur = NaN(1,10000);
    stdA_sur = NaN(1,10000);
    sumA_sur = NaN(1,10000);
    for sur=1:10000
        % clc; disp(sur);
        [MA_sur(1,sur),stdA_sur(1,sur),NAs_surr] = BT_Analysis_NegAsy(Time_Steps_Surr(sur,:),Time_Trigs_Obs);
        sumA_sur(1,sur) = nansum(abs(NAs_surr));
    end
    pmna = sort(abs(MA_sur));
    pmna = find(pmna>=MA_obs,1,'first')/length(pmna);
    psdna = sort(stdA_sur);
    pstdna = find(psdna>=stdA_obs,1,'first')/length(psdna);
    psna = sort([sumA_sur sumNA_obs]);
    psumna = find(psna>=sumNA_obs,1,'first')/length(psna);
    p = [pmna,pstdna,psumna];
    
    sur_order = sumA_sur;
    
elseif strcmp(char(type),'SI')
    % Shannon Entropy and Synchronization Index
    [SI_observed] = BT_Analysis_SyncIndexShannonE(Time_Steps_Obs,Time_Trigs_Obs,10,0);
    SI_sur = NaN(1,10000);
    for sur=1:10000
        %clc; disp(sur)
        SI_sur(1,sur) = BT_Analysis_SyncIndexShannonE(Time_Steps_Surr(sur,:),Time_Trigs_Obs,10,0);
    end
    psi_sur = sort([SI_sur SI_observed]); 
    p=1-find(psi_sur>=SI_observed,1,'first')/length(psi_sur);

    sur_order = SI_sur;
elseif strcmp(char(type),'CIRC')
    % Circular
    [Circ_obs] = BT_Analysis_CircStat(Time_Steps_Obs,Time_Trigs_Obs,'Mean',0); R_obs = Circ_obs(4,3); T_obs = Circ_obs(3,3);
    R_surr = NaN(1,10000);
    T_surr = NaN(1,10000);
    for sur=1:10000
        %clc; disp(sur)
        [Circ_surr] = BT_Analysis_CircStat(Time_Steps_Surr(sur,:),Time_Trigs_Obs,'Mean',0); 
        R_surr(1,sur) = Circ_surr(4,3); T_surr(1,sur) = Circ_surr(3,3);
    end
    pr_sur = sort([R_surr R_obs]);
    pr=1-find(pr_sur>=R_obs,1,'first')/length(pr_sur);
    pt_sur = sort([T_surr T_obs]);
    pt=1-find(pt_sur>=T_obs,1,'first')/length(pt_sur);
    p= [pr pt];
    
    sur_order = R_surr;
    
elseif strcmp(char(type),'CROSSCO')
    [~,~,lagMax_obs,rMax_obs] = BT_Analysis_CrossCor(Time_Steps_Obs,Time_Trigs_Obs,fs,'deviants',8,0);
    rMax_surr   = NaN(1,10000);
    lagMax_surr = NaN(1,10000);
    for sur=1:10000
        clc; disp(sur);
        [~,~,lagMax_surr(1,sur),rMax_surr(1,sur)] = BT_Analysis_CrossCor(Time_Steps_Surr(sur,:),Time_Trigs_Obs,fs,'deviants',8,0);
    end
    prmax_sur = sort(rMax_surr);
    prmax=1-find(prmax_sur>=rMax_obs,1,'first')/length(prmax_sur);
    plagmax_sur = sort(lagMax_surr);
    plagmax=1-find(plagmax_sur>=lagMax_obs,1,'first')/length(plagmax_sur);
    p = [prmax plagmax];
    
    sur_order = rMax_surr;
    
end
if fig==1
    
    Time_Step_Surr_min=Time_Steps_Surr(sur_order==min(sur_order),:);
    Time_Steps_Surr_max=Time_Steps_Surr(sur_order==max(sur_order),:);
    
    
    
    
    AllEvents=sort([Time_Steps_Obs Time_Trigs_Obs Time_Step_Surr_min Time_Steps_Surr_max]);
    time=AllEvents(1)-1:1/fs:AllEvents(end)+1;
    figure;
    subplot(3,2,1)
    BT_Figure_TimeIOI(Time_Steps_Obs,Time_Trigs_Obs,time,fs)
    title('Observed data')
    subplot(3,2,3)
    BT_Figure_TimeIOI(Time_Step_Surr_min,Time_Trigs_Obs,time,fs)
    title('Minimal Sync')
    subplot(3,2,5)
    BT_Figure_TimeIOI(Time_Steps_Surr_max,Time_Trigs_Obs,time,fs)
    title('maximal sync')
    
    if strcmp(char(type),'CIRC')
        [~,rads_obs] = BT_Analysis_CircStat(Time_Steps_Obs,Time_Trigs_Obs,'Mean',0);
        [~,rads_min] = BT_Analysis_CircStat(Time_Step_Surr_min,Time_Trigs_Obs,'Mean',0);
        [~,rads_max] = BT_Analysis_CircStat(Time_Steps_Surr_max,Time_Trigs_Obs,'Mean',0);
        
        subplot(3,2,2)
        StepsRads_TargetTempo = rads_obs(:,3);
        StepsRads_TargetTempo(isnan(StepsRads_TargetTempo))=[];
        circ_plot(StepsRads_TargetTempo,'pretty');
        subplot(3,2,4)
        StepsRads_Surrmin   = rads_min(:,3);
        StepsRads_Surrmin(isnan(StepsRads_Surrmin))=[];
        circ_plot(StepsRads_Surrmin,'pretty');
        subplot(3,2,6)
        StepsRads_Surrmax   = rads_max(:,3);
        StepsRads_Surrmax(isnan(StepsRads_Surrmax))=[];
        circ_plot(StepsRads_Surrmax,'pretty');
    elseif strcmp(char(type),'MNA')
        [~,~,As_obs] = BT_Analysis_NegAsy(Time_Steps_Obs,Time_Trigs_Obs);
        [~,~,As_min] = BT_Analysis_NegAsy(Time_Step_Surr_min,Time_Trigs_Obs);
        [~,~,As_max] = BT_Analysis_NegAsy(Time_Steps_Surr_max,Time_Trigs_Obs);
        
        subplot(3,2,2)
        histogram(As_obs)
        subplot(3,2,4)
        histogram(As_min)
        subplot(3,2,6)
        histogram(As_max)
    end
end

end

