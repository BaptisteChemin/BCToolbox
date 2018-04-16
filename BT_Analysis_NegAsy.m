function [MNA,stdNA,NAs] = BT_Analysis_NegAsy(Time_Steps,Time_Trigs)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
   
    % compute the relative phase
    NAs           = [];
    for s=1:length(Time_Steps)
        Trig_Before_Step        = Time_Trigs(find(Time_Trigs<Time_Steps(s),1,'last'));
        Trig_After_Step         = Time_Trigs(find(Time_Trigs>=Time_Steps(s),1,'first'));
        if ~isempty(Trig_Before_Step) && ~isempty(Trig_After_Step)
            Instant_NA_a        = Time_Steps(s)-Trig_After_Step;
            Instant_NA_b        = Time_Steps(s)-Trig_Before_Step;
            Instant_NA_ab       = [Instant_NA_a Instant_NA_b];
            Instant_NA_i        = abs(Instant_NA_ab)==min(abs(Instant_NA_ab));
            Instant_NA          = Instant_NA_ab(Instant_NA_i);
            NAs   = [NAs Instant_NA]; %#ok<*AGROW>
        else
            NAs   = [NAs NaN];
        end
    end
    MNA = nanmean(NAs);
    stdNA = nanstd(NAs);
end

