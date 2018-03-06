% global steps_t steps_idx
% global SL
% global threshold_slowfall
% global threshold_doublestep
% global keep_second
% global method
% global line time
% global Hs1 Hs2 Hsl Hil His

%% parameters that have to be updated in GetStepLatencies_updateplot
method = 1; % method 1 = threshold; method 2 = filter;
SL = [-2 2]; % levels for the falltime function
threshold_slowfall = .2; % threshold for rejection of artifacts: slow falls
threshold_doublestep = .3; % threshold for rejection of artifacts: repetition of a step
keep_second = 1; % in case of a step repetition, keep the second occurence of the step


%% Processing

if method==1
    % get the indexes and the latencies (s) of identified steps
    [~,steps_idx]=falltime(line,'StateLevels',SL,'Tolerance',5);
    [F,steps_t]=falltime(line,time,'StateLevels',SL,'Tolerance',5);
    % remove steps that are due to "artifacts" of slow fall
    idxf = F > threshold_slowfall;
    steps_t(idxf==1)=[];   
    steps_idx(idxf==1)=[];   
elseif method==2
    % get the indexes and the latencies (s) of identified steps
end

% remove steps that are due to "artifacts" of double steps
D=diff(steps_t);
idxd = D < threshold_doublestep; 
if keep_second==1
    idxd = [idxd 0];
else
    idxd = [0 idxd];
end
steps_t(idxd==1)=[];
steps_idx(idxd==1)=[];

%% Accel and Steps Graph
fig = figure('unit','norm','pos',[.01 .05 .98 .85]); hold on
subplot(2,1,1)
if method==1
    % draw the statelevel lines
    SLdown = NaN(length(time),1);
    SLdown(:,1) = SL(1);
    SLup = NaN(length(time),1);
    SLup(:,1) = SL(2);
    Hs1 = plot(time,SLdown,'g'); hold on; % Handle Statlevel 1
    Hs2 = plot(time,SLup,'r'); % Handel Statelevel 2
elseif metho==2
    % plot the filtered signal
    hold on;
end

% Plot the identified steps and the original signal
axesHandles = findall(fig,'type','axes');
pos_pro=get(axesHandles(1),'position');
pos_pro(2)=0.55; pos_pro(4)=0.4;
pos_pro(1)=0.1186; pos_pro(3)=0.7069;
set(axesHandles(1),'position',pos_pro)
hold on;
steps = NaN(length(time),1);
steps(round(steps_idx)-1) = floor(min(line)); steps(round(steps_idx)) = ceil(max(line));
Hsl = plot(time,steps,'color',[.6 .6 .6]); hold on % Hsl = Handle Steps Line
plot(time,line, 'color',[0 0 0]);
% xlim([0 180])
xlabel('time (s)'); ylabel('amplitude');
title('TimeCourse of LowerLimb Accelerometer')

%% IOI Graph
subplot(2,1,2)
axesHandles = findall(fig,'type','axes');
pos_pro=get(axesHandles(1),'position');
pos_pro(2)=0.28; pos_pro(4)=0.2;
pos_pro(1)=0.1186; pos_pro(3)=0.7069;
set(axesHandles(1),'position',pos_pro)

ISI_t       = diff(steps_t);
fancy_time = time;
fancy_step = NaN(1,length(fancy_time));
for ioi=1:length(ISI_t)
    fancy_step(round(steps_idx(ioi)))=ISI_t(ioi);
end

His = plot(fancy_time,fancy_step,'*'); hold on; % His = Handle Ioi Stars
if exist('fillmissing')==2 %#ok<EXIST>
    fancy_step_line = fillmissing(fancy_step,'linear'); % Hil = Handle Ioi Line
    Hil = plot(fancy_time,fancy_step_line,'-r'); 
end
% xlim([0 180])
ylim([0 5])
xlabel('time (s)'); ylabel('IOI');
title('Inter-Step-Interval along time')



