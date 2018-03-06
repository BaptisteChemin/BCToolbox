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
% delete or add a manually selected step
% if exist('button','var')
%     Cursor = datacursormode(fig);
%     CursorInfo = getCursorInfo(Cursor);
%     if strcmp(button,'Delete Step')
%         [~, removeindex] = min(abs(round(steps_idx)-round(CursorInfo.DataIndex)));
%         steps_t(removeindex(1)) = [];
%         steps_idx(removeindex(1)) = [];
%     elseif strcmp(button,'Add Step')
%         test = sort([steps_t;CursorInfo.Position(1)]);
%         steps_t = sort([steps_t;CursorInfo.Position(1)]);
%         steps_idx = sort([steps_idx;CursorInfo.DataIndex]);
%     end
% end


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
steps(round(steps_idx)-1) = floor(min(line)); steps(round(steps_idx)) = ceil(max(line));
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
