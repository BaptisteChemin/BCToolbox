function [time_keep,IOI_Trigs_nodeviant,IOI_Steps_nodeviant,Time_Trigs_nodeviant,Time_Steps_nodeviant,IOI_Trigs_missing,IOI_Steps_missing,Time_Trigs_missing,Time_Steps_missing] = BT_DataTransform_Time2IEI(Time_Steps,Time_Trigs,time,deviants)

IOI_Trigs               = diff(Time_Trigs);
IOI_Steps               = diff(Time_Steps);
IOI_median              = median([IOI_Trigs IOI_Steps]);

fs                      = 1/mean(diff(time));

if strcmp(char(deviants),'removedeviants')
    Idx_Outliers_Trigs          = isoutlier(IOI_Trigs,'gesd','ThresholdFactor', 0.000000000000001);
    Idx_Outliers_Steps          = isoutlier(IOI_Steps,'gesd','ThresholdFactor', 0.00000000000001);
    
    Time_before_Outliers_Trigs  = Time_Trigs(Idx_Outliers_Trigs);
    Time_before_Outliers_Steps  = Time_Steps(Idx_Outliers_Steps);
    Time_after_Outliers_Trigs   = Time_Trigs([false Idx_Outliers_Trigs(1:end-1)]);
    Time_after_Outliers_Steps   = Time_Steps([false Idx_Outliers_Steps(1:end-1)]);
    
    Time_Outliers_Trigs         = sort([Time_before_Outliers_Trigs Time_after_Outliers_Trigs]);
    Idx_groups_Trigs_start      = find(diff(ismember(Time_Trigs, Time_Outliers_Trigs))==1)+1;
    Idx_groups_Trigs_end        = find(diff(ismember(Time_Trigs, Time_Outliers_Trigs))==-1);
    Idx_groups_Trigs_end        = Idx_groups_Trigs_end+1;
    Time_groups_Trigs_start     = Time_Trigs(Idx_groups_Trigs_start);
    Time_groups_Trigs_end       = Time_Trigs(Idx_groups_Trigs_end);
    [Idx_Trigs_Start,Idx_Trigs_End] = BT_DataTransform_Time2Idx(Time_groups_Trigs_start,Time_groups_Trigs_end,time);
    
    Time_Outliers_Steps         = sort([Time_before_Outliers_Steps Time_after_Outliers_Steps]);
    Idx_groups_Steps_start      = find(diff(ismember(Time_Steps, Time_Outliers_Steps))==1)+1;
    Idx_groups_Steps_end        = find(diff(ismember(Time_Steps, Time_Outliers_Steps))==-1);
    Idx_groups_Steps_end        = Idx_groups_Steps_end+1;
    Time_groups_Steps_start     = Time_Steps(Idx_groups_Steps_start);
    Time_groups_Steps_end       = Time_Steps(Idx_groups_Steps_end);
    [Idx_Steps_Start,Idx_Steps_End] = BT_DataTransform_Time2Idx(Time_groups_Steps_start,Time_groups_Steps_end,time);
    
    time_keep                   = ones(1,length(time));
    IOI_remove                  = NaN(1,length(time));
    Time_remove                 = NaN(1,length(time));
    
    IOI_margin                  = round((IOI_median/4)*fs);
    Time_margin                 = round((IOI_median/2)*fs);
    
    for g=1:length(Idx_Steps_Start)
        time_keep(Idx_Steps_Start(g):Idx_Steps_End(g))= NaN;
        IOI_remove(Idx_Steps_Start(g)+IOI_margin:Idx_Steps_End(g)-IOI_margin) = time(Idx_Steps_Start(g)+IOI_margin:Idx_Steps_End(g)-IOI_margin);
        Time_remove(Idx_Steps_Start(g)+Time_margin:Idx_Steps_End(g)-3*Time_margin) = time(Idx_Steps_Start(g)+Time_margin:Idx_Steps_End(g)-3*Time_margin);
    end
    for g=1:length(Idx_Trigs_Start)
        time_keep(Idx_Trigs_Start(g):Idx_Trigs_End(g))= NaN;
        IOI_remove(Idx_Trigs_Start(g)+IOI_margin:Idx_Trigs_End(g)-IOI_margin) = time(Idx_Trigs_Start(g)+IOI_margin:Idx_Trigs_End(g)-IOI_margin);
        Time_remove(Idx_Trigs_Start(g)+Time_margin:Idx_Trigs_End(g)-3*Time_margin) = time(Idx_Trigs_Start(g)+Time_margin:Idx_Trigs_End(g)-3*Time_margin);
    end
    IOI_remove(isnan(IOI_remove))=[];
    Time_remove(isnan(Time_remove))=[];
    
    
    
    
    IOI_Trigs_nodeviant         = IOI_Trigs;
    IOI_Steps_nodeviant         = IOI_Steps;
    
    Idx_Remove_IOI_Steps = [];
    for s=1:length(Time_Steps)
        if min(abs(IOI_remove-Time_Steps(s)))<.005
            Idx_Remove_IOI_Steps = [Idx_Remove_IOI_Steps s];
        end
    end
    IOI_Steps_nodeviant (Idx_Remove_IOI_Steps-1) = [];
    
    Idx_Remove_IOI_Trigs = [];
    for s=1:length(Time_Trigs)
        if min(abs(IOI_remove-Time_Trigs(s)))<.005
            Idx_Remove_IOI_Trigs = [Idx_Remove_IOI_Trigs s];
        end
    end
    IOI_Trigs_nodeviant (Idx_Remove_IOI_Trigs-1) = [];
    
    
    
    Time_Steps_nodeviant = Time_Steps;
    Idx_Remove_Time_Steps = [];
    for s=1:length(Time_Steps)
        if min(abs(Time_remove-Time_Steps(s)))<.004
            Idx_Remove_Time_Steps = [Idx_Remove_Time_Steps s];
        end
    end
    Time_Steps_nodeviant (Idx_Remove_Time_Steps) = [];
    
    Time_Trigs_nodeviant = Time_Trigs;
    Idx_Remove_Time_Trigs = [];
    for s=1:length(Time_Trigs)
        if min(abs(Time_remove-Time_Trigs(s)))<.004
            Idx_Remove_Time_Trigs = [Idx_Remove_Time_Trigs s];
        end
    end
    Time_Trigs_nodeviant (Idx_Remove_Time_Trigs) = [];
   
    
    % suggest some values to fill missing events
    IOI_Trigs = diff(Time_Trigs);
    IOI_Steps = diff(Time_Steps);
    
    IOI_Steps_corr    = setdiff(IOI_Steps,IOI_Steps_nodeviant);
    IOI_Steps_missing_amount = NaN(size(IOI_Steps_corr));
    for i=1:length(IOI_Steps_corr)
        if IOI_Steps_corr(1,i) > 1.5*IOI_median
            IOI_Steps_missing_amount(1,i) = round(IOI_Steps_corr(1,i)./IOI_median);
        else
            
        end
    end
    
    IOI_Trigs_corr    = setdiff(IOI_Trigs,IOI_Trigs_nodeviant);
    IOI_Trigs_missing_amount =  NaN(size(IOI_Trigs_corr));
    for i=1:length(IOI_Trigs_corr)
        if IOI_Trigs_corr(1,i) > 1.5*IOI_median
            IOI_Trigs_missing_amount(1,i) = round(IOI_Trigs_corr(1,i)./IOI_median);
        else
            
        end
    end
    
    IOI_Steps_missing = min(IOI_Steps_nodeviant) + (max(IOI_Steps_nodeviant)-min(IOI_Steps_nodeviant)).*rand(nansum(IOI_Steps_missing_amount),1);
    IOI_Trigs_missing = min(IOI_Trigs_nodeviant) + (max(IOI_Trigs_nodeviant)-min(IOI_Trigs_nodeviant)).*rand(nansum(IOI_Trigs_missing_amount),1);
    
    Time_Steps_missing = NaN(1,length(IOI_Steps_missing));
    idii=0;
    for id = 1:length(IOI_Steps_corr)
        idx_ioi_dev = IOI_Steps == IOI_Steps_corr(1,id);
        rounds = round(IOI_Steps_corr(1,id)./IOI_median)-1; % because number of step missing == number of IOI missing - 1
        time_start = Time_Steps(idx_ioi_dev);
        for idg = 1:rounds
            idii=idii+1;
            add = min(IOI_Steps_nodeviant) + (max(IOI_Steps_nodeviant)-min(IOI_Steps_nodeviant)).*randn(1,1);
            Time_Steps_missing(1,idii) = time_start + add;
            time_start = time_start + add;
        end
    end
    Time_Steps_missing(isnan(Time_Steps_missing)) = [];
    
    Time_Trigs_missing = NaN(1,length(IOI_Trigs_missing));
    idii=0;
    for id = 1:length(IOI_Trigs_corr)
        idx_ioi_dev = IOI_Trigs == IOI_Trigs_corr(1,id);
        rounds = round(IOI_Trigs_corr(1,id)./IOI_median)-1; % because number of step missing == number of IOI missing - 1
        time_start = Time_Trigs(idx_ioi_dev);
        for idg = 1:rounds
            idii=idii+1;
            add = min(IOI_Trigs_nodeviant) + (max(IOI_Trigs_nodeviant)-min(IOI_Trigs_nodeviant)).*randn(1,1);
            Time_Trigs_missing(1,idii) = time_start + add;
            time_start = time_start + add;
        end
    end
    Time_Trigs_missing(isnan(Time_Trigs_missing)) = [];
    
    if length(IOI_Steps) > length(IOI_Steps_nodeviant) + length(IOI_Steps_missing)
        IOI_Steps_missing = [IOI_Steps_missing setdiff(IOI_Steps,IOI_Steps_nodeviant)];
    end
    
    if length(IOI_Trigs) > length(IOI_Trigs_nodeviant) + length(IOI_Trigs_missing)
        IOI_Trigs_missing = [IOI_Trigs_missing setdiff(IOI_Trigs,IOI_Trigs_nodeviant)];
    end
    
    
    
    
    
    
    
    
    
    
    
else
    time_keep = ones(size(time));
    IOI_Trigs_nodeviant         = IOI_Trigs;
    IOI_Steps_nodeviant         = IOI_Steps;
end