function [R,rads] = BT_Analysis_CircStat(Time_Steps,Time_Trigs,SelfMethod,fig)

%XPGAIT_CircStat performs circular statistics on steps latencies and
%targets latencies.
%   R is a 4x3 matrix.
%   Lines:
%   (1) is IOI used for computation of Rads_Steps
%   (2) is Reyglig Test (nonhomogeneous distribution if <.05);
%   (3) is tetha angle in radian;
%   (4) is result vector R.
%   Columns:
%   (1) is data aligned to mean self tempo
%   (2) is data aligned to mean stimulus tempo
%   (3) is data aligned to instantaneous stimulus tempo
% A PRIORI, only columns (1) and (3) are needed 

%   rads is a numberofsteps matrix X 3
%   Columns:
%   (1) is data aligned to mean self tempo
%   (2) is data aligned to mean stimulus tempo
%   (3) is data aligned to instantaneous stimulus tempo

if fig == 1
    figure;
    set(gcf,'Units','normalized')
    pos = get(gcf,'position');
    pos(1,3)=1; pos(1,1)=0;
    set(gcf,'position',pos)
end
%% Dispertion around self mean tempo
if strcmp(char(SelfMethod),'Best')
    IOI_Steps_central = BT_Analysis_CircStat_adjust(Time_Steps,0);
    Rads_Steps_SelfMean = 2*pi*(Time_Steps/IOI_Steps_central);
elseif strcmp(char(SelfMethod),'Mean')
    IEI = diff(Time_Steps);
    Idx_Outliers_IEI          = isoutlier(IEI,'gesd');
    IEI(Idx_Outliers_IEI)     = [];
    Rads_Steps_SelfMean = 2*pi*(Time_Steps/mean(IEI));
elseif strcmp(char(SelfMethod),'Median') 
    Rads_Steps_SelfMean = 2*pi*(Time_Steps/median(diff(Time_Steps)));
end
rads(:,1) = Rads_Steps_SelfMean';
R(1,1) = mean(diff(Time_Steps));
R(2,1) = circ_rtest(Rads_Steps_SelfMean');
R(3,1) = circ_mean(Rads_Steps_SelfMean');
R(4,1) = circ_r(Rads_Steps_SelfMean');
if fig == 1
    subplot(1,3,1)
    circ_plot(Rads_Steps_SelfMean','pretty');
    title('Dispertion around mean self tempo')
    text(.7,.9,['R = ' num2str(R(4,1))])
end

if ~isempty (Time_Trigs)
    %% Dispertion around mean auditory target
    Rads_Steps_TriggerMean = 2*pi*(Time_Steps/mean(diff(Time_Trigs)));
    rads(:,2) = Rads_Steps_TriggerMean';
    R(1,2) = mean(diff(Time_Trigs));
    R(2,2) = circ_rtest(Rads_Steps_TriggerMean');
    R(3,2) = circ_mean(Rads_Steps_TriggerMean');
    R(4,2) = circ_r(Rads_Steps_TriggerMean');
    if fig == 1
        subplot(1,3,2)
        circ_plot(Rads_Steps_TriggerMean','pretty');
        title('Dispertion around mean auditory target')
        text(.7,.9,['R = ' num2str(R(4,2))])
    end
    %% Dispertion around real auditory target
    Rads_Steps_TriggerInstant   = [];
    for s=1:length(Time_Steps)
        Trig_Before_Step        = Time_Trigs(find(Time_Trigs<Time_Steps(s),1,'last'));
        Trig_After_Step         = Time_Trigs(find(Time_Trigs>=Time_Steps(s),1,'first'));
        if ~isempty(Trig_Before_Step) && ~isempty(Trig_After_Step)
            Instant_Period      = Trig_After_Step-Trig_Before_Step;
            Rads_Steps_TriggerInstant   = [Rads_Steps_TriggerInstant 2*pi*((Time_Steps(s)-Trig_Before_Step)/Instant_Period)]; %#ok<*AGROW>
        else
            Rads_Steps_TriggerInstant   = [Rads_Steps_TriggerInstant NaN];
        end
    end
    rads(:,3) = Rads_Steps_TriggerInstant';
    Rads_Steps_TriggerInstant(isnan(Rads_Steps_TriggerInstant)) = [];
    R(1,3) = mean(diff(Time_Trigs));
    R(2,3) = circ_rtest(Rads_Steps_TriggerInstant');
    R(3,3) = circ_mean(Rads_Steps_TriggerInstant');
    R(4,3) = circ_r(Rads_Steps_TriggerInstant');
    if fig == 1
        subplot(1,3,3)
        circ_plot(Rads_Steps_TriggerInstant','pretty');
        title('Dispertion around real auditory target')
        text(.7,.9,['R = ' num2str(R(4,3))])
    end
else
    rads(:,2) = NaN(length(rads),1);
    
    R(1,2) = NaN;
    R(2,2) = NaN;
    R(3,2) = NaN;
    R(4,2) = NaN;
    
    rads(:,3) = NaN(length(rads),1);
    
    R(1,3) = NaN;
    R(2,3) = NaN;
    R(3,3) = NaN;
    R(4,3) = NaN;    
end

%%%%%%%%%

%%%%%%%%%%%%%% OLD VERSION WITH WARPING AND OTHER STUFFS;
%function [R,rads,IOI_Steps_central,IOI_Steps_central_warp] = BT_Analysis_CircStat(Time_Steps,Time_Trigs,fs,warpisneeded)
%   R is a 4x4 matrix.
%   Lines:
%   (1) is IOI used for computation of Rads_Steps
%   (2) is Reyglig Test (nonhomogeneous distribution if <.05);
%   (3) is tetha angle in radian;
%   (4) is result vector R.
%   Columns:
%   (1) is non warped data aligned to mean self tempo
%   (2) is non warped data aligned to mean stimulus tempo
%   (3) is warped data aligned to mean self tempo
%   (4) is warped data aligned to instantaneous stimulus tempo
% A PRIORI, only columns (1) and (4) are needed (case the stimulation is non-stationnary), and columns (1) and (2) case the stimulation is stationary are needed.
% But we perform all the possibilities of analysis "just in case".
%   rads is a numberofsteps matrix X 4
%   Columns:
%   (1) is non warped data aligned to mean self tempo
%   (2) is non warped data aligned to mean stimulus tempo
%   (3) is warped data aligned to mean self tempo
%   (4) is warped data aligned to instantaneous stimulus tempo
% Time_Steps(Time_Steps<=0) = [];
% figure;
% set(gcf,'Units','normalized')
% pos = get(gcf,'position');
% pos(1,3)=1; pos(1,1)=0;
% set(gcf,'position',pos)

% % Dispertion around self mean tempo
% % FIND BEST IOI_STEPS
% IOI_Steps_central = BT_Analysis_CircStat_adjust(Time_Steps);
% % Voila
% Rads_StepsSelf = 2*pi*(Time_Steps/IOI_Steps_central);
% rads(:,1) = Rads_StepsSelf';
% subplot(1,4,1)
% circ_plot(Rads_StepsSelf','pretty');
% title('Dispertion around self tempo')
% R(1,1) = IOI_Steps_central;
% R(2,1) = circ_rtest(Rads_StepsSelf');
% R(3,1) = circ_mean(Rads_StepsSelf');
% R(4,1) = circ_r(Rads_StepsSelf');
% 
% % Dispertion around the auditory target
% IOI_Trigs = mean(diff(Time_Trigs));
% Rads_Steps = 2*pi*(Time_Steps/IOI_Trigs);
% rads(:,2) = Rads_Steps';
% subplot(1,4,2)
% circ_plot(Rads_Steps','pretty');
% title('Dispertion around mean auditory target')
% R(1,2) = IOI_Trigs;
% R(2,2) = circ_rtest(Rads_Steps');
% R(3,2) = circ_mean(Rads_Steps');
% R(4,2) = circ_r(Rads_Steps');
% 
% if warpisneeded == 1
%     time            = (Time_Trigs(1)-1):1/fs:(Time_Trigs(end)+1);
%     Line_Trigs      = zeros(size(time));
%     Line_Steps      = zeros(size(time));
%     
%     [idx_steps,idx_trigs] = BT_DataTransform_Time2Idx(Time_Steps,Time_Trigs,time);
%     
%     Line_Trigs(idx_trigs)   = 1; Line_Trigs(idx_trigs+1) = 1; Line_Trigs(idx_trigs+2) = 1; Line_Trigs(idx_trigs+3) = 1; Line_Trigs(idx_trigs+4) = 1;
%     Line_Steps(idx_steps)   = 1; Line_Steps(idx_steps+1) = 1; Line_Steps(idx_steps+2) = 1; Line_Steps(idx_steps+3) = 1; Line_Steps(idx_steps+4) = 1;
%     Line_Trigs              = Line_Trigs(1:length(time));
%     Line_Steps              = Line_Steps(1:length(time));
%     
%     t_samples     = Time_Trigs;
%     t_query       = t_samples(1):mean(diff(Time_Trigs)):t_samples(end)+1;
%     t_query       = t_query(1:length(t_samples));
%     
%     [Line_Trigs_warp]       = BT_PreProcess_TimeWarping(time,t_samples,t_query,Line_Trigs);
%     [Line_Steps_warp]       = BT_PreProcess_TimeWarping(time,t_samples,t_query,Line_Steps);
%     
%     Line_Trigs_warp(Line_Trigs_warp~=0) = 1;
%     Line_Steps_warp(Line_Steps_warp~=0) = 1;
%     
%     [~,Time_Trigs_warp] = risetime(Line_Trigs_warp,time);
%     Time_Trigs_warp = Time_Trigs_warp';
%     [~,Time_Steps_warp] = risetime(Line_Steps_warp,time);
%     Time_Steps_warp = Time_Steps_warp';
%     
%     % Dispertion around self mean tempo
%     IOI_Steps_central_warp = BT_Analysis_CircStat_adjust(Time_Steps_warp);
%     Rads_StepsSelf_warp = 2*pi*(Time_Steps_warp/IOI_Steps_central_warp);
%     rads(:,3) = Rads_StepsSelf_warp';
%     subplot(1,4,3)
%     circ_plot(Rads_StepsSelf_warp','pretty');
%     title('Dispertion around warped self tempo')
%     R(1,3) = IOI_Steps_central_warp;
%     R(2,3) = circ_rtest(Rads_StepsSelf_warp');
%     R(3,3) = circ_mean(Rads_StepsSelf_warp');
%     R(4,3) = circ_r(Rads_StepsSelf_warp');
%     
%     % Dispertion around the auditory target
%     IOI_Trigs_warp = mean(diff(Time_Trigs_warp));
%     Rads_Steps_warp = 2*pi*(Time_Steps_warp/IOI_Trigs_warp);
%     rads(:,4) = Rads_Steps_warp';
%     subplot(1,4,4)
%     circ_plot(Rads_Steps_warp','pretty');
%     title('Dispertion around real auditory target')
%     R(1,4) = IOI_Trigs_warp;
%     R(2,4) = circ_rtest(Rads_Steps_warp');
%     R(3,4) = circ_mean(Rads_Steps_warp');
%     R(4,4) = circ_r(Rads_Steps_warp');
%     
% else
%     R(:,3:4) = NaN;
%     rads(:,3:4) = NaN;
%     %IOI_Trigs_warp = NaN;
% end
%%%%%%%%%%%%%% OLD VERSION WITH WARPING AND OTHER STUFFS;

end

