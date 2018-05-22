% global steps_t steps_idx
% global line time
% global Hsl Hil His

%% Processing
% delete or add a manually selected step
if exist('button','var')
    Cursor = datacursormode(fig);
    CursorInfo = getCursorInfo(Cursor);
    if strcmp(button,'Delete Step')
        [~, removeindex] = min(abs(round(steps_idx)-round(CursorInfo.DataIndex)));
        steps_t(removeindex(1)) = [];
        steps_idx(removeindex(1)) = [];
    elseif strcmp(button,'Add Step')
        %test = sort([steps_t CursorInfo.Position(1)]);
        steps_t = sort([steps_t CursorInfo.Position(1)]);
        steps_idx = sort([steps_idx CursorInfo.DataIndex]);
    end
end

%% Accel and Steps Graph

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
