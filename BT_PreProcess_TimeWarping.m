function [out_data] = BT_PreProcess_TimeWarping(fs,index_samples,index_query,data)
% 
% 
%
%
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


data_length     = length(data);
t               = (1/fs:1/fs:1/fs*data_length);

out_data=zeros(size(data));


if strcmp(mode,'interm') % this try is not making any difference. I need to think about the area of influence of the frontier points...
    index_samples_interm  = NaN(1,length(index_samples)-1);
    for idx = 1:length(index_samples_interm)
        index_samples_interm(idx) = mean(index_samples(idx:idx+1));
    end
    index_query_interm = NaN(1,length(index_query)-1);
    for idx = 1:length(index_query_interm)
        index_query_interm(idx) = mean(index_query(idx:idx+1));
    end
    x=[t(1),index_query_interm(2:end),t(end)];
    y=[t(1),index_samples_interm(2:end),t(end)];

elseif strcmp(mode,'straight')

    x=[t(1),index_query(2:end),t(end)];
    y=[t(1),index_samples(2:end),t(end)];

end

t1=interp1(x,y,t);
out_data=interp1(t,data',t1)';
