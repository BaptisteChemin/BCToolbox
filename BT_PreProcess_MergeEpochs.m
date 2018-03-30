ids = {'01 EEG2';'02';'03 EEG2';'04 EEG2';'05 EEG2';'06 EEG2';'07 EEG2';'08';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20'};
for pat=1:length(ids)
    type        = 'Patient';   % 'Patient'; 'Control'; 'Test';
    id          = char(ids(pat));%'01'
    load(['ep XPPARKMOTPRIM ' char(type(1)) num2str(id) '.mat']); % = data
    load(['ep XPPARKMOTPRIM ' char(type(1)) num2str(id) '.lw6'], '-mat'); % = header
    out_data=data; out_header=header;
    %process data
    clc
    disp(['*** Processing : ' char(type(1)) num2str(id(1:2))]);
    for epoch_idx = 1:header.datasize(1)
        load(['epoch ' num2str(epoch_idx) ' XPPARKMOTPRIM ' char(type(1)) num2str(id) '.mat']); % = data
    
        out_data(epoch_idx,:,:,:,:,:) = data;
        clear data;
    end
    data = out_data;
    header = out_header;
    save (['epoch XPPARKMOTPRIM ' char(type(1)) num2str(id) '.mat'], 'data')
    save (['epoch XPPARKMOTPRIM ' char(type(1)) num2str(id) '.lw6'], 'header')
end