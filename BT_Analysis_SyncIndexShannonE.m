function [SI] = BT_Analysis_SyncIndexShannonE(Time_Steps,Time_Trigs,BinSize,fig)

%XPGAIT_SyncIndexShannonE returns a Synchronization Index based on a measure of Shannon Entropy in the relative phase between phase of stim and phase of movement.
% SI = 1-SE/ln N'
% The synchronization index ranges from 0, where the spreading of relative
% phase is maximal (all phases lie in different bins), to 1, when a delta
% function like probability distribution is found (all phases lie in a
% single bin).

%% Comfortable CircStat visualisation
% figure;
% Rads_Steps = 2*pi*(Time_Steps/mean(diff(Time_Trigs)));
% subplot(1,2,1)
% circ_plot(Rads_Steps','pretty');
% title('Dispertion around mean auditory target')
% text(.7,.9,['R = ' num2str(circ_r(Rads_Steps'))])
% 
% Rads_Steps_Instant   = [];
% for s=1:length(Time_Steps)
%     Trig_Before_Step        = Time_Trigs(find(Time_Trigs<Time_Steps(s),1,'last'));
%     Trig_After_Step         = Time_Trigs(find(Time_Trigs>=Time_Steps(s),1,'first'));
%     if ~isempty(Trig_Before_Step) && ~isempty(Trig_After_Step)
%         Instant_Period      = Trig_After_Step-Trig_Before_Step;
%         Rads_Steps_Instant   = [Rads_Steps_Instant 2*pi*((Time_Steps(s)-Trig_Before_Step)/Instant_Period)]; %#ok<*AGROW>
%     end
% end
% subplot(1,2,2)
% circ_plot(Rads_Steps_Instant','pretty');
% title('Dispertion around real auditory target')
% text(.7,.9,['R = ' num2str(circ_r(Rads_Steps_Instant'))])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% ?rel(n)=?music(n)-?motion(n)
% relative phase, in degrees.
Deg_Rel_TrigSteps           = [];
for s=1:length(Time_Steps)
    Trig_Before_Step        = Time_Trigs(find(Time_Trigs<Time_Steps(s),1,'last'));
    Trig_After_Step         = Time_Trigs(find(Time_Trigs>=Time_Steps(s),1,'first'));
    if ~isempty(Trig_Before_Step) && ~isempty(Trig_After_Step)
        Instant_Period      = Trig_After_Step-Trig_Before_Step;
        Instant_Phase       = circ_rad2ang(2*pi*((Time_Steps(s)-Trig_Before_Step)/Instant_Period));
        if 180 <= Instant_Phase && Instant_Phase < 360 % keep instant phase within a -180 to 180 ° range
            Instant_Phase = Instant_Phase-360;
        end
        Deg_Rel_TrigSteps   = [Deg_Rel_TrigSteps Instant_Phase]; %#ok<*AGROW>
        %Rad_Rel_TrigSteps   = [Deg_Rel_TrigSteps 2*pi*((Time_Steps(s)-Trig_Before_Step)/Instant_Period)]; %#ok<*AGROW>
    else 
        Deg_Rel_TrigSteps   = [Deg_Rel_TrigSteps NaN];
    end
end
%% Figure representing the data and the relative phase
if fig == 1
    time = (Time_Trigs(1)-1):1/250:(Time_Trigs(end)+1);
    figure;
    subplot(2,1,1)
    BT_Figure_TimeIOI(Time_Steps,Time_Trigs,time,250)
    subplot(2,1,2)
    plot(Deg_Rel_TrigSteps); ylim([-180 180]);   % Adjust it to the real time scale, not bin number.
end

Deg_Rel_TrigSteps(isnan(Deg_Rel_TrigSteps)) = []; % because I added a NaN for the figure.

%% Compute the Shannon Entropy of the data

xbins = -180:BinSize:180;
[n,g]=hist(Deg_Rel_TrigSteps,xbins);%create bin repartition for a bin size of 10
%[n,g]=hist(Deg_Rel_TrigSteps,BinSize);
for i=1:length(n)
    p(i)=n(i)/sum(n);% sum(n) == length(Deg_Rel_TrigSteps);  %associate a probability for the bin
end
p=p(p~=0); %keep non zero elements
ent=-sum(p.*log(p));

%% Compute the Synchronization Index by comparing SE of the data to maximal possible entropy 

nmax_indiv_values = sum(n)/length(n);
nmax = nmax_indiv_values.*ones(size(n));
for i=1:length(nmax)
    pmax(i)=nmax(i)/sum(nmax);%associate a probability for the bin
end
pmax=pmax(pmax~=0); %keep non zero elements
entmax=-sum(pmax.*log(pmax));

SI = 1-(ent/entmax);
end

