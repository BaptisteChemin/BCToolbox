function [Signal_Steps,Signal_Trigs,time] = BT_DataTransform_Time2Signal(Time_Steps,Time_Trigs,fs,type)
    All_Events              = sort([Time_Steps Time_Trigs]);
    time                    = (All_Events(1)-1):1/fs:(All_Events(end)+1);
    
    Signal_Steps            = zeros(size(time));
    Signal_Trigs            = zeros(size(time));
    [idx_steps,idx_trigs]   = BT_DataTransform_Time2Idx(Time_Steps,Time_Trigs,time);
    if strcmp(type,'trigger_type')
        Signal_Steps(idx_steps) = 1; Signal_Steps(idx_steps+1) = 1; Signal_Steps(idx_steps+2) = 1; Signal_Steps(idx_steps+3) = 1; Signal_Steps(idx_steps+4) = 1; 
        Signal_Trigs(idx_trigs) = 1; Signal_Trigs(idx_trigs+1) = 1; Signal_Trigs(idx_trigs+2) = 1; Signal_Trigs(idx_trigs+3) = 1; Signal_Trigs(idx_trigs+4) = 1;
    elseif strcmp(type,'wave_type')
        for i=1:length(idx_steps)-1
            wave_d  = idx_steps(i+1)-idx_steps(i);
            wave_t  = 0:wave_d;
            wave_f  = 1/(wave_d);
            wave_y  = -1*cos(2*pi*wave_f*wave_t);
            Signal_Steps(idx_steps(i):idx_steps(i+1))=wave_y+1;
        end
        for i=1:length(idx_trigs)-1
            wave_d  = idx_trigs(i+1)-idx_trigs(i);
            wave_t  = 0:wave_d;
            wave_f  = 1/(wave_d);
            wave_y  = -1*cos(2*pi*wave_f*wave_t);
            Signal_Trigs(idx_trigs(i):idx_trigs(i+1))=wave_y+1;
        end
    end
end

