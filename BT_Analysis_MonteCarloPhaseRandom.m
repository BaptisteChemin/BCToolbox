function [Time_Steps_Surr] = BT_Analysis_MonteCarloPhaseRandom(Time_Steps,Time_Trigs,fs,nsurr,type)
%XPGAIT_MonteCarloPhaseRandom generates surrogate Data from time Serie (Prichard 1994), and statistically test phase synchronisation of actual data by comparing it to 10000 surrogates.
% adapted from (c) 2011, Carlos
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

%% pre-processing
% transform time series into usefull signal (serie of Inter-Event-Intervals or
% succession of waveforms)

All_Events              = sort([Time_Trigs Time_Steps]);
medianperiod            = median([diff(Time_Steps) diff(Time_Trigs)]);
time                    = (All_Events(1)-medianperiod/2):1/fs:(All_Events(end)+medianperiod/2);

% Deal with missing events (part 1)
% sometimes, some events were missing, generating artifacually large
% IOI. Here, I compensate so there are no such artifacts (that would
% disturb the phase randomisation!). Later, the times of events will
% be shifted to adjust acording to this compensation.
[time_keep,IOI_Trigs_nodeviant,IOI_Steps_nodeviant,~,~,IOI_Trigs_missing,IOI_Steps_missing,Time_Trigs_missing,Time_Steps_missing] = BT_DataTransform_Time2IEI(Time_Steps,Time_Trigs,time,'removedeviants');

IOI_Trigs = diff(Time_Trigs);
IOI_Steps = diff(Time_Steps);
IOI_Steps_corr    = setdiff(IOI_Steps,IOI_Steps_nodeviant);
IOI_Trigs_corr    = setdiff(IOI_Trigs,IOI_Trigs_nodeviant);

% generate the "Signals" as either sine waves or series of Inter Event
% Intervals. Fills missing events.
if strcmp(char(type),'wave_type')
    Time_Steps_nomissing = sort([Time_Steps Time_Steps_missing]);
    Time_Trigs_nomissing = sort([Time_Trigs Time_Trigs_missing]);
    [Signal_Steps,Signal_Trigs,time] = BT_DataTransform_Time2Signal(Time_Steps_nomissing,Time_Trigs_nomissing,fs,'wave_type');
    Signal_Steps                = Signal_Steps';
    Signal_Trigs                = Signal_Trigs';
elseif strcmp(char(type),'IOI_type')
    Signal_Steps                = [IOI_Steps_nodeviant';IOI_Steps_missing];
    Signal_Trigs                = [IOI_Trigs_nodeviant';IOI_Trigs_missing];
end

% crop the data so it is suitable for spectral manipulations
[nfrms,nts]                 = size(Signal_Steps);
if rem(nfrms,2) == 0
    nfrms                   = nfrms-1;
    Signal_Steps            = Signal_Steps(1:nfrms,:);
    if strcmp(char(type),'wave_type')
        time                    = time(:,1:nfrms);
    elseif strcmp(char(type),'IOI_type')
        last_step = mean(Signal_Steps) + (max(Signal_Steps)-min(Signal_Steps))/2*randn;
    end
end 

% Get parameters
len_ser                     = (nfrms-1)/2;
interv1                     = 2:len_ser+1;
interv2                     = len_ser+2:nfrms;

%% Phase randomisation
% Fourier transform of the original dataset
Signal_Steps_fft            = fft(Signal_Steps);
    
% Create the surrogate data one by one
Signal_Steps_Surr           = zeros(nfrms, nsurr);
Time_Steps_Surr             = zeros(nsurr, length(Time_Steps));
    
    for k = 1:nsurr
        ph_rnd                   = rand([len_ser 1]);   %%%%  RANDOM NUMBER BETWEEN 0 AND 1
        
        % Create the random phases
        ph_interv1               = repmat(exp( 2*pi*1i*ph_rnd),1,nts);
        ph_interv2               = conj( flipud( ph_interv1));
        
        % Randomize the time serie
        Signal_Steps_fft_surr             = Signal_Steps_fft;
        Signal_Steps_fft_surr(interv1,:)  = Signal_Steps_fft(interv1,:).*ph_interv1;
        Signal_Steps_fft_surr(interv2,:)  = Signal_Steps_fft(interv2,:).*ph_interv2;
        
        % Inverse transform
        Signal_Steps_Surr(:,k)   = real(ifft(Signal_Steps_fft_surr));
        
        %% post-processing
        % convert the randomized signal into data of the original kind.
        if strcmp(char(type),'wave_type')
            Signal_Steps_temp           = -Signal_Steps_Surr(:,k);
            Signal_Steps_temp(isnan(time_keep)) = min(Signal_Steps_temp);
            [~,Time_Steps_Surr_temp]    = findpeaks(Signal_Steps_temp,time);
            if length(Time_Steps_Surr_temp)>=length(Time_Steps)
                Time_Steps_Surr(k,:)     = Time_Steps_Surr_temp(1:length(Time_Steps));
            elseif length(Time_Steps_Surr_temp)<length(Time_Steps)
                Time_Steps_Surr(k,:)     = [Time_Steps_Surr_temp NaN(1,length(Time_Steps)-length(Time_Steps_Surr_temp))];
            end
            
        elseif strcmp(char(type),'IOI_type')
            first_step = Time_Steps(1) + mean(Signal_Steps_Surr(:,k))/2*randn;
            
            % shift Time_Steps in the case there was some artefactually large IEI.
            for id = 1:length(IOI_Steps_missing)
                idx_dev = find(IOI_Steps == IOI_Steps_corr(id,1));
                var = Signal_Steps_Surr(idx_dev,k)-mean(Signal_Steps_Surr(:,k));
                Signal_Steps_Surr(idx_dev,k) = IOI_Steps_corr(id,1)+var;
            end

            if exist('last_step','var')
                Time_Steps_Surr(k,:)     = cumsum([first_step; Signal_Steps_Surr(:,k);last_step]);
            else
                Time_Steps_Surr(k,:)     = cumsum([first_step; Signal_Steps_Surr(:,k)]);
            end            
        end
        
          
    end   
end

