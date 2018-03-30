function [Time_Steps_Surr] = BT_Analysis_MonteCarloPhaseRandom(Time_Steps,Time_Trigs,fs,nsurr)
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


[Signal_Steps,Signal_Trigs,time] = BT_DataTransform_Time2Signal(Time_Steps,Time_Trigs,fs,'wave_type');
Signal_Steps                = Signal_Steps';
Signal_Trigs                = Signal_Trigs';

[nfrms,nts]                 = size(Signal_Steps);
if rem(nfrms,2) == 0
    nfrms                   = nfrms-1;
    Signal_Steps            = Signal_Steps(1:nfrms,:);
    time                    = time(:,1:nfrms);
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
   ph_rnd                   = rand([len_ser 1]);
   
   % Create the random phases for all the time series
   ph_interv1               = repmat(exp( 2*pi*1i*ph_rnd),1,nts);
   ph_interv2               = conj( flipud( ph_interv1));
   
   % Randomize all the time series simultaneously
   IOI_fft_surr             = Signal_Steps_fft;
   IOI_fft_surr(interv1,:)  = Signal_Steps_fft(interv1,:).*ph_interv1;
   IOI_fft_surr(interv2,:)  = Signal_Steps_fft(interv2,:).*ph_interv2;
   
   % Inverse transform
   Signal_Steps_Surr(:,k)   = real(ifft(IOI_fft_surr));  
   [~,Time_Steps_Surr_temp] = findpeaks(-Signal_Steps_Surr(:,k),time);
   Time_Steps_Surr(k,:)     = Time_Steps_Surr_temp(1:length(Time_Steps));
end

