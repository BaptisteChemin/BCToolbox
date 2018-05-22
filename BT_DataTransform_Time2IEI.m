function [time_keep,IOI_Trigs_nodeviant,IOI_Steps_nodeviant] = BT_DataTransform_Time2IEI(Time_Steps,Time_Trigs,time,deviants)

IOI_Trigs               = diff(Time_Trigs);
IOI_Steps               = diff(Time_Steps);
IOI_median              = median([IOI_Trigs IOI_Steps]);

fs                      = 1/mean(diff(time));

if strcmp(char(deviants),'removedeviants')
    Idx_Outliers_Trigs          = isoutlier(IOI_Trigs,'gesd');
    Idx_Outliers_Steps          = isoutlier(IOI_Steps,'gesd');
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
    Time_remove                 = NaN(1,length(time));
    
    margin                      = round((IOI_median/4)*fs);
    for g=1:length(Idx_Steps_Start)
        time_keep(Idx_Steps_Start(g):Idx_Steps_End(g))= NaN;
        Time_remove(Idx_Steps_Start(g)+margin:Idx_Steps_End(g)-margin) = time(Idx_Steps_Start(g)+margin:Idx_Steps_End(g)-margin);
    end
    for g=1:length(Idx_Trigs_Start)
        time_keep(Idx_Trigs_Start(g):Idx_Trigs_End(g))= NaN;
        Time_remove(Idx_Steps_Start(g)+margin:Idx_Steps_End(g)-margin) = time(Idx_Steps_Start(g)+margin:Idx_Steps_End(g)-margin);
    end
    
    
    IOI_Trigs_nodeviant         = IOI_Trigs;
    IOI_Steps_nodeviant         = IOI_Steps;
    
    Time_remove(isnan(Time_remove))=[];
    Time_Steps_nodeviant = Time_Steps;
    Idx_Remove_Steps = [];
    for s=1:length(Time_Steps_nodeviant)
        if min(abs(Time_remove-Time_Steps_nodeviant(s)))<.004
            Idx_Remove_Steps = [Idx_Remove_Steps s];
        end
    end
    %Time_Steps_nodeviant(Idx_Remove_Steps) = [];
    IOI_Steps_nodeviant (Idx_Remove_Steps-1) = [];
    
    
    Time_Trigs_nodeviant = Time_Trigs;
    Idx_Remove_Trigs = [];
    for s=1:length(Time_Trigs_nodeviant)
        if min(abs(Time_remove-Time_Trigs_nodeviant(s)))<.004
            Idx_Remove_Trigs = [Idx_Remove_Trigs s];
        end
    end
    %Time_Trigs_nodeviant(Idx_Remove_Trigs) = [];
    IOI_Trigs_nodeviant (Idx_Remove_Trigs-1) = [];
    
else
    time_keep = time;
    IOI_Trigs_nodeviant         = IOI_Trigs;
    IOI_Steps_nodeviant         = IOI_Steps;
end