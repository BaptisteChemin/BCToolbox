function [Times_AC, Times_NC, Times_I, Envelope_AC, Envelope_NC, Envelope_I, Shuffle_O] = BT_GenSignal_ACNC(numberofevents,meanperiod,jitmax,jitstd,genaudio)
%% preamble
fs                  = 44100;
frequax             = 1000/meanperiod/2*(2/numberofevents:2/numberofevents:2);
if mod(numberofevents,2) == 0
    numberofevents = numberofevents+1;
    disp('The number of events has to be odd, therefore I added 1 to your original NoE...')
end
% defaultfilenames    =  {'Isochronous';
%     'NonAutocorrelated';
%     'Autocorrelated'};
%% generate the acoustic event
slope_t         = (meanperiod-meanperiod*(jitmax+0.02))/1000;
slope_ts        = (0:1/fs:slope_t-1/fs);
noteenv         = zeros(1,length(slope_ts));
for i=1:length(slope_ts)
    t           = (i)/(length(slope_ts)/2);
    noteenv(i)  = exp(-t*3)*(1-exp(-t*20))*2;
end
noteenv         = noteenv/(max(noteenv)+0.01);
slope_on        = noteenv(1:find(noteenv==max(noteenv)));
slope_off       = noteenv(find(noteenv==max(noteenv)):end);

carry1=rand(1,5*size(noteenv,2));  carry2=carry1.^2; carry3=carry2-1; carry4=carry3.^2;
carry5=2*carry4; carry=carry5-1; clear carry1; clear carry2; clear carry3; clear carry4; clear carry5;% create a uniform white noise

%% generate the signals
%%%%%%%% ISOCHRONOUS %%%%%%%%
Isochronous_IOI_s           = NaN(1,numberofevents); Isochronous_IOI_s(1,:) = meanperiod/1000;
Isochronous_envelope        = [];
Isochronous_index_time      = (1/fs:meanperiod/1000:numberofevents*meanperiod/1000); 
Isochronous_index_samples   = round(Isochronous_index_time*fs);
Isochronous_sequence        = [];
clicktrig_I                 = [];

    jitter = (Isochronous_IOI_s(1) - slope_ts(length(slope_on)) - slope_ts(length(slope_off)));
    jitter_add = NaN(1,round(jitter*fs)); 
    jitter_add(1:end) = slope_off(end);
    env = [slope_on, slope_off, jitter_add];
    BIP = carry(1:length(env)).*env;
    
for bip=1:length(Isochronous_IOI_s)
    Isochronous_sequence(1,round(Isochronous_index_samples(1,bip)):(round(Isochronous_index_samples(1,bip))+length(BIP)-1))=BIP;
    Isochronous_envelope(1,round(Isochronous_index_samples(1,bip)):(round(Isochronous_index_samples(1,bip))+length(BIP)-1))=env;
    clicktrig_I(1,round(Isochronous_index_samples(1,bip)):(round(Isochronous_index_samples(1,bip))+0.01*fs)-1)=ones(0.01*fs,1);
end
clicktrig_I = [clicktrig_I zeros(1,length(Isochronous_sequence)-length(clicktrig_I))];


% generate the jitters for auto-correlated and NonAutocorrelated sequences
distribution_check = 1;
while distribution_check == 1
    jitter_raw          = randn(numberofevents-1,1);
    jitter_mz           = jitter_raw-mean(jitter_raw);
    
    jitter_f            = fft(jitter_mz);
    mask                = ones(length(jitter_f)/2,1);
    mask                = mask./frequax(1:length(mask))';
    mask                = [mask; zeros(length(mask)-1,1)];
    mask(ceil(length(mask)/2)+1:end) = flip(mask(1:floor(length(mask)/2)));
    jitter_ff           = jitter_f;
    jitter_ff(2:end)    = jitter_ff(2:end).*mask;

    jitter_ac           = real(ifft(jitter_ff));
    jitter_n            = jitter_ac/max(abs(jitter_ac));
    
    std_check = 1;
    jitter_pc           = jitter_n;
    while std_check == 1
        var = std(jitter_pc);
        if var >= jitstd-0.0005 && var <= jitstd+0.0005
            std_check = 0;
        elseif var<jitstd
            jitter_pc           = jitter_pc*1.001;
        elseif var>jitstd
            jitter_pc           = jitter_pc*0.999;
        end
    end    
    if jitmax < max(abs(jitter_pc))
        indi = find(abs(jitter_pc)>jitmax);
        for ji=1:length(indi)
            if jitter_pc(indi(ji))>0
                jitter_pc(indi(ji))    = jitmax;
            elseif  jitter_pc(indi(ji))<0
                jitter_pc(indi(ji))    = -jitmax;
            end
        end
    end  
    
    IOI_s   = ((jitter_pc+1)*meanperiod/1000)';
    distribution_check              = +(adtest(IOI_s)); clc
    if distribution_check == 1
        close all;
    end
end

Autocorrelated_jitter_pc        = jitter_pc;
Autocorrelated_IOI_s            = ((Autocorrelated_jitter_pc+1)*meanperiod/1000)';
NonAutocorrelated_jitter_pc     = Autocorrelated_jitter_pc(randperm(numberofevents-1));
NonAutocorrelated_IOI_s         = ((NonAutocorrelated_jitter_pc+1)*meanperiod/1000)';

% Shufle_O can be used to shuflle a vector corresponding to IOI AC to make it correspond to IOI NC
[sorted_AC_IOI,order_ac_1]                  = sort(Autocorrelated_IOI_s);
[~,order_ac]                    = sort(order_ac_1);
[sorted_NC_IOI,order_nc_1]                  = sort(NonAutocorrelated_IOI_s);
[~,order_nc]                    = sort(order_nc_1);
Shuffle_O = 1:length(order_ac); 
Shuffle_O = Shuffle_O(order_ac_1);
Shuffle_O = Shuffle_O(order_nc);
%%%%%%%%%% Autocorrelated %%%%%%%%
Autocorrelated_envelope         = [];
Autocorrelated_index_time       = [1/fs cumsum(Autocorrelated_IOI_s)];          % here we add the first beat (on first sample, 1/fs); to keep!
Autocorrelated_index_samples    = Autocorrelated_index_time*fs;
Autocorrelated_sequence         = [];
clicktrig_AC                    = [];
for bip = 1:length(Autocorrelated_IOI_s)
    jitter                      = (Autocorrelated_IOI_s(bip) - slope_ts(length(slope_on)) - slope_ts(length(slope_off)));
    jitter_add                  = NaN(1,round(jitter*fs)); 
    jitter_add(1:end)           = slope_off(end);
    
    env                         = [slope_on, slope_off, jitter_add];
    BIP                         = carry(1:length(env)).*env;
    Autocorrelated_sequence(1,round(Autocorrelated_index_samples(1,bip)):(round(Autocorrelated_index_samples(1,bip))+length(BIP)-1))=BIP;
    Autocorrelated_envelope(1,round(Autocorrelated_index_samples(1,bip)):(round(Autocorrelated_index_samples(1,bip))+length(BIP)-1))=env;
    clicktrig_AC(1,round(Autocorrelated_index_samples(1,bip)):(round(Autocorrelated_index_samples(1,bip))+0.01*fs)-1)=ones(0.01*fs,1);
end
clicktrig_AC = [clicktrig_AC zeros(1,length(Autocorrelated_sequence)-length(clicktrig_AC))];

%%%%%%%% NonAutocorrelated %%%%%%%%
NonAutocorrelated_envelope      = [];
NonAutocorrelated_index_time    = [1/fs cumsum(NonAutocorrelated_IOI_s)];          % here we add the first beat (on first sample, 1/fs); to keep!
NonAutocorrelated_index_samples = NonAutocorrelated_index_time*fs;
NonAutocorrelated_sequence      = [];
clicktrig_NC                    = [];
for bip = 1:length(NonAutocorrelated_IOI_s)
    %clc; disp(bip);
    jitter                      = (NonAutocorrelated_IOI_s(bip) - slope_ts(length(slope_on)) - slope_ts(length(slope_off)));
    jitter_add                  = NaN(1,round(jitter*fs)); 
    jitter_add(1:end)           = slope_off(end);
    
    env                         = [slope_on, slope_off, jitter_add];
    BIP                         = carry(1:length(env)).*env;
    NonAutocorrelated_sequence(1,round(NonAutocorrelated_index_samples(1,bip)):(round(NonAutocorrelated_index_samples(1,bip))+length(BIP)-1))=BIP;
    NonAutocorrelated_envelope(1,round(NonAutocorrelated_index_samples(1,bip)):(round(NonAutocorrelated_index_samples(1,bip))+length(BIP)-1))=env;
    clicktrig_NC(1,round(NonAutocorrelated_index_samples(1,bip)):(round(NonAutocorrelated_index_samples(1,bip))+0.01*fs)-1)=ones(0.01*fs,1);
end
clicktrig_NC                    = [clicktrig_NC zeros(1,length(NonAutocorrelated_sequence)-length(clicktrig_NC))];

blank   = zeros(1,5*fs);
blank2  = zeros(1,15*fs);

%% SAVE
Times_AC    = Autocorrelated_index_time;
Times_NC    = NonAutocorrelated_index_time;
Times_I     = Isochronous_index_time;
Envelope_AC = Autocorrelated_envelope;
Envelope_NC = NonAutocorrelated_envelope;
Envelope_I  = Isochronous_envelope;

if genaudio == 1
    audiowrite('Autocorrelated.wav', [[blank Autocorrelated_sequence blank2]' [blank clicktrig_AC blank2]'],fs);
    audiowrite('NonAutocorrelated.wav', [[blank NonAutocorrelated_sequence blank2]' [blank clicktrig_NC blank2]'],fs);
    audiowrite('Isochronous.wav', [[blank Isochronous_sequence blank2]' [blank clicktrig_I blank2]'],fs);
end