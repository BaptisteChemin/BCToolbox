function [MNA,stdNA,As,Ad] = BT_Analysis_NegAsy(Times_Step,Times_Trig)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
   

    As                          = []; % Asynchronies in second
    Ad                          = []; % Asynchronies in degree
    
    
    for s=1:length(Times_Step)
        Trig_Before_Step        = Times_Trig(find(Times_Trig<Times_Step(s),1,'last'));
        Trig_After_Step         = Times_Trig(find(Times_Trig>=Times_Step(s),1,'first'));
        if ~isempty(Trig_Before_Step) && ~isempty(Trig_After_Step)
            %% Asynchrony in seconds
            Instant_NA_a        = Times_Step(s)-Trig_After_Step;
            Instant_NA_b        = Times_Step(s)-Trig_Before_Step;
            Instant_NA_ab       = [Instant_NA_a Instant_NA_b];
            Instant_NA_i        = abs(Instant_NA_ab)==min(abs(Instant_NA_ab)); % select the value closest to 
            Instant_NA          = Instant_NA_ab(Instant_NA_i);
            As                  = [As Instant_NA]; %#ok<*AGROW>
            
            %% Asynchrony in degree
            Instant_Period      = Trig_After_Step-Trig_Before_Step;
            Instant_Phase       = circ_rad2ang(2*pi*((Times_Step(s)-Trig_Before_Step)/Instant_Period));
            if 180 <= Instant_Phase && Instant_Phase < 360 % keep instant phase within a -180 to 180 ° range
                Instant_Phase = Instant_Phase-360;
            end
            Ad   = [Ad Instant_Phase]; %#ok<*AGROW>
            
        else
            As   = [As NaN];
            Ad   = [Ad NaN];
        end
    end
    MNA = nanmean(As);
    stdNA = nanstd(As);
end

