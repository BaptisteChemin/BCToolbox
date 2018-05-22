part = {'P02';'P04';'P08';'P09';'P10';'P11';'P12';'P14';'P15';'P16';'P17';'P18';'P19';'P20'}; 
F = .5;
frequi = F:F:6*F;

for p = 1:length(part)
    load(['TOPO OK Filtered XPPARKMOTPRIM ' char(part(p)) '.mat']); % = data
    load(['TOPO OK Filtered XPPARKMOTPRIM ' char(part(p)) '.lw6'], '-mat');
    
    % size(data)

    Fb = round((frequi/header.xstep)+1);
    
    chunk=zeros(size(data,2),1);
    for ch =1:size(data,2)
        line = squeeze(data(1,ch,1,1,1,:));
        chunk(ch) = sum(line(Fb));
    end
    
    data(1,:,1,1,1,end)=chunk;
    % plot(squeeze(data(1,:,1,1,1,end)))
    
    save (['TOPO OK Filtered XPPARKMOTPRIM ' char(part(p)) '.mat'], 'data')
    save (['TOPO OK Filtered XPPARKMOTPRIM ' char(part(p)) '.lw6'], 'header')


end