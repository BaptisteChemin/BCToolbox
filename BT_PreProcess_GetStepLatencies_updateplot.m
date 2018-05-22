% global steps_t steps_idx
% global SL sl1 sl2
% global threshold_slowfall sf
% global threshold_doublestep sr
% global keep_second srp
% global method
% global line time
% global Hs1 Hs2 Hsl Hil His

%% get the updated variables (if updated)
if exist('sl1','var')
    SL(1) = sl1;
end
if exist('sl2','var')
    SL(2) = sl2;
end
if exist('sf','var')
    threshold_slowfall = sf;
end
if exist('sr','var')
    threshold_doublestep = sr;
end
if exist('srp','var')
    keep_second = srp;
end
if exist('ar','var')
    adapt_range = ar;
end


%% Processing
if method==1
    [F,steps_t]=falltime(line,time,'StateLevels',SL,'Tolerance',5);
    [~,steps_idx]=falltime(line,'StateLevels',SL,'Tolerance',5);
    % delete the slow fall artefactual steps
    idxf = F > threshold_slowfall;
    steps_t(idxf==1)=[];   
    steps_idx(idxf==1)=[];
end
% delete the repetition artefactual steps
D=diff(steps_t);
idxd = D < threshold_doublestep; 
if keep_second==1
    idxd = [idxd 0];
else
    idxd = [0 idxd];
end
steps_t(idxd==1)=[];
steps_idx(idxd==1)=[];
steps_idx = round(steps_idx);


% adjust the steps to the minimal value within a_r range of bins
if adapt_range > 0
    for step=1:length(steps_idx)
        if length(line)>steps_idx(step)+adapt_range && steps_idx(step)-adapt_range > 0
            [negpeakvalue,negpeaktime]=findpeaks(-line(steps_idx(step)-adapt_range:steps_idx(step)+adapt_range),time(steps_idx(step)-adapt_range:steps_idx(step)+adapt_range));
            negpeaktime = negpeaktime(find(negpeakvalue==max(negpeakvalue)));
        else
            negpeaktime = [];
        end
        if ~isempty(negpeaktime)
            steps_t(step)   = negpeaktime;
            steps_idx(step) = find(time==negpeaktime);
        end
    end
end


%% Accel and Steps Graph
if method==1
    SLdown = NaN(length(time),1);
    SLdown(:,1) = SL(1);
    SLup = NaN(length(time),1);
    SLup(:,1) = SL(2);
    set(Hs1,'ydata',SLdown)
    set(Hs2,'ydata',SLup)
end

steps = NaN(length(time),1);
steps(round(steps_idx)-1) = floor(min(line)-.2); steps(round(steps_idx)) = ceil(max(line)+.2);
set(Hsl,'ydata',steps)

%% IOI Graph
ISI_t       = diff(steps_t);
fancy_step = NaN(1,length(fancy_time));
for ioi=1:length(ISI_t)
    fancy_step(round(steps_idx(ioi)))=ISI_t(ioi);
end
if exist('fillmissing')==2 %#ok<EXIST>
    fancy_step_line = fillmissing(fancy_step,'linear'); % Hil = Handle Ioi Line
    set(Hil,'ydata',fancy_step_line)
end
set(His,'ydata',fancy_step)
