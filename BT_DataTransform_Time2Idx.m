function [idx_steps,idx_trigs] = BT_DataTransform_Time2Idx(Time_Steps,Time_Trigs,time)
    %time            = (Time_Trigs(1)-1):1/fs:(Time_Trigs(end)+1);
    All_Events = sort([Time_Steps Time_Trigs]);
    if ~isempty(All_Events(All_Events<=time(1)))
        error('Some time values are outside the range of time vector.')
    end
    
    idx_trigs       = NaN(size(Time_Trigs));
    idx_steps       = NaN(size(Time_Steps));
    
    for t = 1:length(Time_Trigs)
        idx_trigs(t)= find(abs(time - Time_Trigs(t))==min(abs(time - Time_Trigs(t))),1,'first');
    end
    for t = 1:length(Time_Steps)
        idx_steps(t)= find(abs(time - Time_Steps(t))==min(abs(time - Time_Steps(t))),1,'first');
    end
end

