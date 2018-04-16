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

% transform time series into usefull signal
if strcmp(char(type),'wave_type')
    [Signal_Steps,Signal_Trigs,time] = BT_DataTransform_Time2Signal(Time_Steps,Time_Trigs,fs,'wave_type');
    Signal_Steps                = Signal_Steps';
    Signal_Trigs                = Signal_Trigs';
elseif strcmp(char(type),'IOI_type')
    Signal_Steps                = diff(Time_Steps)';
    Signal_Trigs                = diff(Time_Trigs)';
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
    
% Fourier transform of the original dataset
Signal_Steps_fft            = fft(Signal_Steps);
    
% Create the surrogate recording blocks one by one
Signal_Steps_Surr           = zeros(nfrms, nsurr);
Time_Steps_Surr             = zeros(nsurr, length(Time_Steps));
    
    for k = 1:nsurr
        %clc; disp(k);
        ph_rnd                   = rand([len_ser 1]);   %%%%  RANDOM NUMBER BETWEEN 0 AND 1
        %ph_rnd                   = randn([len_ser 1]);
        
        % Create the random phases for all the time series
        ph_interv1               = repmat(exp( 2*pi*1i*ph_rnd),1,nts);
        ph_interv2               = conj( flipud( ph_interv1));
        
        % Randomize all the time series simultaneously
        Signal_Steps_fft_surr             = Signal_Steps_fft;
        Signal_Steps_fft_surr(interv1,:)  = Signal_Steps_fft(interv1,:).*ph_interv1;
        Signal_Steps_fft_surr(interv2,:)  = Signal_Steps_fft(interv2,:).*ph_interv2;
        
        % Inverse transform
        Signal_Steps_Surr(:,k)   = real(ifft(Signal_Steps_fft_surr));
        
        
        if strcmp(char(type),'wave_type')
            [~,Time_Steps_Surr_temp] = findpeaks(-Signal_Steps_Surr(:,k),time);
            if length(Time_Steps_Surr_temp)>=length(Time_Steps)
                Time_Steps_Surr(k,:)     = Time_Steps_Surr_temp(1:length(Time_Steps));
            elseif length(Time_Steps_Surr_temp)<length(Time_Steps)
                Time_Steps_Surr_temp     = [Time_Steps_Surr_temp NaN(1,length(Time_Steps)-length(Time_Steps_Surr_temp))];
                Time_Steps_Surr(k,:)     = Time_Steps_Surr_temp;
            end
        elseif strcmp(char(type),'IOI_type')
            first_step = Time_Steps(1) + mean(Signal_Steps_Surr(:,k))/2*randn;
            if exist('last_step','var')
                Time_Steps_Surr(k,:)     = cumsum([first_step; Signal_Steps_Surr(:,k);last_step]);
            else
                Time_Steps_Surr(k,:)     = cumsum([first_step; Signal_Steps_Surr(:,k)]);
            end
        end
    end   
end

