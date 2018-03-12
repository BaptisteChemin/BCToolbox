function [out_data] = BT_PreProcess_TimeWarping(time,time_samples,time_query,data)

% Author : 
% Baptiste Chemin
% Institute of Neurosciences (IONS)
% Universite catholique de louvain (UCL)
% Belgium
% 
% Contact : baptiste.chemin@uclouvain.be
% This subfunction is part of WarpModeling project
%
mode = 'straight';

% if length(time)>1
%     time = time;
% elseif length(time) == 1
%     data_length     = length(data);
%     fs=time;
%     time=1/fs:1/fs:data_length/fs;
% end

if strcmp(mode,'interm') % this try is not making any difference. I need to think about the area of influence of the frontier points...
    index_samples_interm  = NaN(1,length(time_samples)-1);
    for idx = 1:length(index_samples_interm)
        index_samples_interm(idx) = mean(time_samples(idx:idx+1));
    end
    index_query_interm = NaN(1,length(time_query)-1);
    for idx = 1:length(index_query_interm)
        index_query_interm(idx) = mean(time_query(idx:idx+1));
    end
    x=[time(1),index_query_interm(2:end),time(end)];
    y=[time(1),index_samples_interm(2:end),time(end)];

elseif strcmp(mode,'straight')

    x=[time(1),time_query(2:end),time(end)];
    y=[time(1),time_samples(2:end),time(end)];

end

t1=interp1(x,y,time);
out_data=interp1(time,data',t1)';
