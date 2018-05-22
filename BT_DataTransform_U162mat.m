%% CONVERSION AND ADJUSTMENTS FROM.U16 to .mat

clear;
% load the .wav and .U16 files
[wavname, wavpath] = uigetfile('*.wav', 'Get the stimulation .wav file');
clc; disp (wavname); 
eval(sprintf('%s',['audio = audioread(''' char(wavpath) char(wavname) ''');']));

[u16name, U16path] = uigetfile('*.U16', 'Get the raw .U16 file');
curdir = cd;
eval(sprintf('%s',['cd ''' char(U16path) ''''])); % needs to be in the proper directory to open the U16 file

fid = fopen(u16name,'r');
if fid < 0
    error(['Cannot open file "' u16name '".']);
end
n = floor(length(fread(fid)) / 2);
fseek(fid,0,'bof');
gait = fread(fid,n,'uint16')';
fclose(fid);

eval(sprintf('%s',['cd ''' char(curdir) ''''])); % needs to be in the proper directory to open the U16 file

% management of buffer overflow
flowing         = [0 gait(1:11:end)];
flowing_N       = flowing/max(abs(flowing));
flowing_D       = [NaN diff(flowing(2:end)) NaN];
target          = 1:length(flowing);

flowing_continuous = flowing;
flowing_continuous(1:2)=1:2;

[a,b]=find(flowing_D==-65535);

for r=1:length(b)
    flowing_continuous(b(r)+1:end) = flowing_continuous(b(r)+1:end)+abs(-65535-1);
end
flowing_continuous(3:end)=flowing_continuous(3:end)+2;


reorder = flowing_continuous;
reorder(flowing_continuous)=reorder;

delete = find(reorder~=target);



% separate the channels
trigger         = [0 gait(11:11:end)];
trigger_N       = trigger/max(abs(trigger));
trigger_N(trigger_N>=.25)=1;
trigger_N(trigger_N<.25)=0;
trigger_N(2:20) =1;

sound           = [0 gait(2:11:end)];
jambe_gauche_x  = [0 gait(3:11:end)]; jambe_gauche_x  = jambe_gauche_x - mean(jambe_gauche_x);
jambe_gauche_y  = [0 gait(4:11:end)]; jambe_gauche_y  = jambe_gauche_y - mean(jambe_gauche_y);
jambe_droite_x  = [0 gait(7:11:end)]; jambe_droite_x  = jambe_droite_x - mean(jambe_droite_x);
jambe_droite_y  = [0 gait(8:11:end)]; jambe_droite_y  = jambe_droite_y - mean(jambe_droite_y);
main_gauche_x   = [0 gait(5:11:end)]; main_gauche_x   = jambe_gauche_x - mean(jambe_gauche_x);
main_gauche_y   = [0 gait(6:11:end)]; main_gauche_y   = jambe_gauche_y - mean(jambe_gauche_y);
main_droite_x   = [0 gait(9:11:end)]; main_droite_x   = jambe_droite_x - mean(jambe_droite_x);
main_droite_y   = [0 gait(10:11:end)]; main_droite_y   = jambe_droite_y - mean(jambe_droite_y);
clear gait; clear filename; clear fid; clear n;
    
% correct or suppress signal when buffer overflow occurs
trigger_N(flowing_continuous)=trigger_N;
trigger_N(delete) = 0;
jambe_gauche_x(flowing_continuous)=jambe_gauche_x;
jambe_gauche_x(delete) = NaN;
jambe_gauche_y(flowing_continuous)=jambe_gauche_y;
jambe_gauche_y(delete) = NaN;
jambe_droite_x(flowing_continuous)=jambe_droite_x;
jambe_droite_x(delete) = NaN;
jambe_droite_y(flowing_continuous)=jambe_droite_y;
jambe_droite_y(delete) = NaN;

% crop the segment of interest according to start and stop indices
figure;
hold on
plot(trigger_N); plot(flowing_N);
flowing_D(flowing_D==-65535) =1;
plot(flowing_D);
hold off
title('Set the segment of interest limits (start and stop indices):')
legend('Triggers','Floow','DF')
ylim([0 1.2]);
start = str2double(input('start index:  ','s'));
stop  = str2double(input('stop index:  ','s'));
close
jambe_gauche_x  = jambe_gauche_x(start:stop);
jambe_gauche_y  = jambe_gauche_y(start:stop);
jambe_droite_x  = jambe_droite_x(start:stop);
jambe_droite_y  = jambe_droite_y(start:stop);
main_gauche_x   = main_gauche_x(start:stop);
main_gauche_y   = main_gauche_y(start:stop);
main_droite_x   = main_droite_x(start:stop);
main_droite_y   = main_droite_y(start:stop);
sound           = sound(start:stop);
trigger_N       = trigger_N(start:stop);


figure;
plot(audio(:,2));
title('Set the segment of interest limits (start and stop indices):')
legend('Triggers audio')
ylim([0 1.2]);
start = str2double(input('start index:  ','s'));
stop  = str2double(input('stop index:  ','s'));
close
audio  = audio(start:stop,:);



% 
[lvl]=risetime([0 trigger_N]);
recordednumberoftrig = length(lvl); clear lvl;

[~,start_sample, stop_sample] = risetime(audio(:,2));
start = floor(start_sample(1)); stop = ceil(stop_sample(end)); clear start_sample; clear stop_sample;

audio_crop = resample(audio(start-20:stop+20,2)',250,44100); audio_crop = [0 audio_crop max(audio_crop) 0]; clear start; clear stop;
audio_crop(audio_crop>=.25)=1;
audio_crop(audio_crop<.25)=0;
[lvl] = risetime(audio_crop);
realnumberoftrig = length(lvl); clear lvl;


if realnumberoftrig ~= recordednumberoftrig %|| realnumberoftrig < restart % no problem of not-recorded trigger
    warning('Uncoherent number of performed vs recorded triggers');
    
    if ~isempty(delete)
        error('You need to continue manually and adjust the selection zone in audio signal. Beware of buffer overflow problems in recorded data.')
    else
        [~,start_sample, stop_sample] = risetime(audio(:,2),'statelevels',[0.005 0.05]);
        start = floor(start_sample(1)); stop = ceil(stop_sample(recordednumberoftrig)); clear start_sample; clear stop_sample;
        
        audio_crop = resample(audio(start-20:stop+20,2)',250,44100); audio_crop = [0 audio_crop max(audio_crop) 0]; clear start; clear stop;
        audio_crop(audio_crop>=.25)=1;
        audio_crop(audio_crop<.25)=0;
        [lvl] = risetime(audio_crop);
        realnumberoftrig = length(lvl); clear lvl;
    end
end



[~,start_sample, stop_sample]=risetime([0 trigger_N]);
start_sample = start_sample-1; if floor(start_sample(1))==0; start_sample(1) =1; end
stop_sample  = stop_sample-1;

start = floor(start_sample(1)); stop = ceil(stop_sample(end)); clear start_sample; clear stop_sample;
trigger_crop = trigger_N(start:stop+1);
jambe_gauche_x_crop = jambe_gauche_x(start:stop+1);
jambe_gauche_y_crop = jambe_gauche_y(start:stop+1);
jambe_droite_x_crop = jambe_droite_x(start:stop+1);
jambe_droite_y_crop = jambe_droite_y(start:stop+1);


% adjust the size of signals

lac = length(audio_crop);
ltc = length(trigger_crop);

ratio_pre = ltc/lac;
factor = floor(sqrt(2^31/ratio_pre));
ratio = floor(ratio_pre*factor);

trigger = resample(trigger_crop,factor,ratio);
jambe_gauche_x = resample(jambe_gauche_x_crop,factor,ratio);
jambe_gauche_y = resample(jambe_gauche_y_crop,factor,ratio);
jambe_droite_x = resample(jambe_droite_x_crop,factor,ratio);
jambe_droite_y = resample(jambe_droite_y_crop,factor,ratio);
 
%%PLOT IT%%
close all
audio_goodsize = [audio_crop NaN NaN NaN];
audio_goodsize = audio_goodsize(1:length(trigger));
trigger(trigger>=.25)=1;
trigger(trigger<.25)=0;

figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
hold on;
plot(trigger);
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
plot(audio_goodsize*max(trigger));
hold off;
legend('recorded_trigger', 'real_trigger')
subplot(2,1,2)
hold on;
plot(trigger);
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
plot(audio_goodsize*max(trigger));
hold off;
xlim([length(trigger)-1500 length(trigger)])


check = questdlg('Are the recordings well alligned with the real stims?','Quality Check','yes','no','no');
% save if everything is fine
if strcmp(char(check),'yes')
    % eval(sprintf('%s',['filename=''' char(type(1)) num2str(id) '_GAIT_BEHAVIORAL_' char(stimuli(s)) '_' char(date) '.mat'';']));
    uisave ({'trigger', 'jambe_droite_x', 'jambe_droite_y', 'jambe_gauche_x', 'jambe_gauche_y'},char(u16name(1:end-4)));
end
close all

