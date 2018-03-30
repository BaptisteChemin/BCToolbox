ids = {'01 EEG2';'02';'03 EEG2';'04 EEG2 lastpremissing';'05 EEG2';'06 EEG2';'07 EEG2';'08';'09';'10';'11';'12';'14';'15';'16';'17';'18';'19';'20'};
for pat=1:length(ids)
    type        = 'Patient';   % 'Patient'; 'Control'; 'Test';
    id          = char(ids(pat));%'01'
    
    load(['dc pre-post ep XPPARKMOTPRIM ' char(type(1)) num2str(id) '.mat']); % = data
    load(['dc pre-post ep XPPARKMOTPRIM ' char(type(1)) num2str(id) '.lw6'], '-mat'); % = header
    
    %process data
    clc
    disp(['*** Processing : ' char(type(1)) num2str(id(1:2))]);
    in_header   = header;
    in_data     = data;
    clear header; clear data;
    
    for epoch_idx = 1:size(in_data,1)
        [header,data]=RLW_arrange_epochs(in_header,in_data,epoch_idx);
        save (['epoch ' num2str(epoch_idx) ' dc pre-post XPPARKMOTPRIM ' char(type(1)) num2str(id) '.mat'], 'data')
        save (['epoch ' num2str(epoch_idx) ' dc pre-post XPPARKMOTPRIM ' char(type(1)) num2str(id) '.lw6'], 'header')
        clear header; clear data;
    end
end