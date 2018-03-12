function [central_IOI] = BT_Analysis_CircStat_adjust(Time_Onsets)

IOI_min     = min(diff(Time_Onsets));
IOI_max     = max(diff(Time_Onsets));
res         = 0.001; 

IOI_test    = IOI_min:res:IOI_max;
R           = NaN(size(IOI_test));

for t = 1:length(IOI_test)
    Rads_StepsSelf = 2*pi*(Time_Onsets/IOI_test(t));
    R(1,t) = circ_r(Rads_StepsSelf');
end
central_IOI = IOI_test(1,R==max(R));