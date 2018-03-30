% function [out_header,out_data] = BT_PreProcess_Remove_OffsetDeviants(header,data,varargin)
%%%% This function interpolates outliers values in order to removes offsets artifacts. 
%%% IT NEEDS MATLAB 2016b OR HIGHER VERSION OF MATLAB.

% ex: [BT_PreProcess_Remove_OffsetDeviants]=LW_warping(header,data,'window',window,'threshold',threshold,'bins',bins);
ids = {'01 EEG2';'02';'03 EEG2';'04 EEG2 lastpremissing';'05 EEG2';'06 EEG2';'07 EEG2';'08';'09';'10';'11';'12';'14';'15';'16';'17';'18';'19';'20'};
for pat=1:length(ids)
    type        = 'Patient';   % 'Patient'; 'Control'; 'Test';
    id          = char(ids(pat));%'01';        % number of the participant
    
    load(['pre-post ep XPPARKMOTPRIM ' char(type(1)) num2str(id) '.mat']); % = data
    load(['pre-post ep XPPARKMOTPRIM ' char(type(1)) num2str(id) '.lw6'], '-mat'); % = header
    
    %process data
    %clc
    disp(['*** Processing : ' char(type(1)) num2str(id(1:2))]);
    fs=1/header.xstep;
    %parse varargin
    if  ~exist('varargin','var') || isempty(varargin)
        window = round(.5 * fs);
        thresh = 10;
        bins = 2;
    else
        %window
        a=find(strcmpi(varargin,'window'));
        if isempty(a)
            window = round(.5 * fs);
        else
            window = varargin{a+1};
        end
        %threshold
        a=find(strcmpi(varargin,'threshold'));
        if isempty(a)
            thresh = 10;
        else
            thresh = varargin{a+1};
        end
        %bins
        a=find(strcmpi(varargin,'bins'));
        if isempty(a)
            bins = 2;
        else
            bins = varargin{a+1};
        end
    end
    
    %%
    out_data = data;  
    for t=1:size(data,1)
        idxs_real = [];
        for ch=1:size(data,2)
            line = double(squeeze(out_data(t,ch,1,1,1,:)));
            for b=1:window:length(line)-window
                win = line(b:b+window);
                win_detrend = detrend(win);
                win_nomax = win_detrend; win_nomax(win_nomax==max(abs(win_nomax)))=NaN; win_nomax(win_nomax==max(abs(win_nomax)))=NaN; % I removed some of the max values to compute std, so we have a better estimate even tough there might be some super deviant values in the set.
                m=nanmedian(win_nomax);
                s=nanstd(win_nomax);
                clear win_nomax;
                % note that the 2017b version of matlab introduced a new
                % function "isoutlier" that might be very helpfull here.
                th = [m-thresh*s m+thresh*s];
                idxs_n=find(win_detrend<th(1));
                idxs_p=find(win_detrend>th(2));
                idxs = sort([idxs_n; idxs_p])';
                if ~isempty (idxs)
                    idxs_real = [idxs_real (idxs+b-1)];
                end
                clear idxs; clear idxs_n; clear idxs_p;
            end
            %         end
        end
        idxs_real=unique(sort(idxs_real));
        mult = [1 diff(idxs_real)];
        multiples=[1 find(mult~=1) length(mult)+1];
        clear mult;
        %%%%%%%
        % Like this, it will interpolate the signal in all EEG channels, even
        % if the artifact was detected only on a few channels. I did it because
        % it happens that the artifact spread out to other channels, even tough
        % the amplitude is then quite limitated. To be discussed!
        %%%%%%%
        for ch=1:size(data,2)
            line = double(squeeze(data(t,ch,1,1,1,:)));
            if ~isempty(idxs_real)
                disp(num2str(ch))
                for g = 1:length(multiples)-1
                    idx = idxs_real(multiples(g):multiples(g+1)-1);
                    idx = idx(1)-bins:idx(end)+bins;
                    line(idx) = NaN;
                end                
            end
           xax = 1/fs:1/fs:length(line)/fs;
           line_corr = fillmissing(line,'linear','SamplePoints',xax);
           out_data(t,ch,1,1,1,:) = single(line_corr);
        end
        clear multiples;
    end
    
    data=out_data;
    
    header.name = [header.name, 'JumpFilt'];
    save (['JumpFilt pre-post XPPARKFLUCT ' char(type(1)) num2str(id) '.mat'], 'data')
    save (['JumpFilt pre-post XPPARKFLUCT ' char(type(1)) num2str(id) '.lw6'], 'header')
    
end