part = {'P04';'P08';'P09';'P10';'P11';'P12';'P14';'P15';'P16';'P17';'P18';'P19';'P20'}; 
cond = {'Pre';'Post'};

for p = 1:length(part)
    load(['tapping XPPARKMOTPRIM ' char(part(p)) '.mat']); % = data
    load(['tapping XPPARKMOTPRIM ' char(part(p)) '.lw6'], '-mat'); % = header
    fs = 1/header.xstep;
    for c = 1:length(cond)
        line = double(squeeze(data(c,1,1,1,1,:)));
        plot(line)
        direction = input('what is the direction of tap?   (up,do):         ','s');
        %direction = questdlg('what is the direction of maximal acceleration?','Channel Orientation','falling','rising','falling');
        if strcmp(direction,'do')
            line = -line;
        end
        time = header.xstart:1/fs:(length(line)/fs)+header.xstart; time=time(1:end-1);
        BT_PreProcess_GetStepLatencies
        jambe_droite_times = steps_t;
    end
end