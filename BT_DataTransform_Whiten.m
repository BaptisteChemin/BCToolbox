function [Signal_Leader_Filt,Signal_Follower_Filt,Signal_Leader_Error,Signal_Follower_Error] = BT_DataTransform_Whiten(Signal_Leader,Signal_Follower,Na)
    
    All_Signals             = sort([Signal_Leader, Signal_Follower]);
    
    Signal_Leader           = (Signal_Leader - mean(All_Signals))/std(All_Signals);
    Signal_Follower         = (Signal_Follower - mean(All_Signals))/std(All_Signals);
    
    
    if ~exist('Na','var')
        Na = 100;
    end
    a                       = lpc(Signal_Leader,Na);
    Signal_Leader_Filt      = fftfilt(a,Signal_Leader);
    % Signal_Leader_Filt = Signal_Leader_Filt(Na+1:end);   
    Signal_Follower_Filt    = fftfilt(a,Signal_Follower);
    
    
    Signal_Leader_Error     = Signal_Leader_Filt - Signal_Leader;
    Signal_Follower_Error   = Signal_Follower_Filt - Signal_Follower;
end

