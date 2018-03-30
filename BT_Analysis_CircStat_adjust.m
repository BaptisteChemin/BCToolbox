function [central_IOI] = BT_Analysis_CircStat_adjust(Time_Onsets,fig)

IOI_min     = min(diff(Time_Onsets));
IOI_max     = max(diff(Time_Onsets));
res         = 0.001; 
fs          = 1/res;
time        = min(Time_Onsets)-1:res:max(Time_Onsets)+1;


IOI_test    = IOI_min:res:IOI_max;
R           = NaN(size(IOI_test));
if fig == 1
    figure;
    h1 = subplot(2,1,1);
    h2 = subplot(2,1,2);
end
for t = 1:length(IOI_test)
    Rads_StepsSelf = 2*pi*(Time_Onsets/IOI_test(t));
    R(1,t) = circ_r(Rads_StepsSelf');
    Time_Trigs = min(Time_Onsets):IOI_test(t):max(Time_Onsets);
    if fig == 1;
        cla(h1)
        h1 = subplot(2,1,1);
        BT_Figure_TimeIOI(Time_Onsets,Time_Trigs,time,fs)
        h2 = subplot(2,1,2);
        circ_plot(Rads_StepsSelf','pretty');
        text(.7,.9,['R = ' num2str(R(1,t))])
        text(.7,1,['P = ' num2str(IOI_test(t))])
        if floor(IOI_test(t)*fs) == floor(mean(diff(Time_Onsets))*fs);
            pause(1)
            text(.7,8,'MEAN PERIOD')
        else
            pause(.05)
        end
    end
end
close
central_IOI = IOI_test(1,R==max(R));