function [step_times_crop,trigger_times_crop,N] = BT_DataTransform_LatenciesDataCropAuto(step_times,trigger_times,fs,start,numofsteps)
%CROP the data: 

%% time
AllEvents = sort([trigger_times step_times]);
time = AllEvents(1)-1:1/fs:AllEvents(end)+1;

start_time          = step_times(start)-.1;%str2double(input('Start Time:  ','s'));
if length(step_times)>=start+numofsteps
    stop_time           = step_times(start+numofsteps)+.1;%str2double(input('Stop Time:  ','s'));
else
    stop_time           = step_times(end);
end

start_step          = find(step_times - start_time>=0,1,'first');
stop_step           = find(step_times - stop_time<=0,1,'last');

start_trig          = find(trigger_times - start_time>=0,1,'first');
stop_trig           = find(trigger_times - stop_time<=0,1,'last');


step_times_crop     = step_times(start_step:stop_step);
if start_trig<stop_trig
    trigger_times_crop  = trigger_times(start_trig:stop_trig);
else
    trigger_times_crop  = [];
end

N(1,1)=length(trigger_times_crop);
N(2,1)=length(step_times_crop);

end

