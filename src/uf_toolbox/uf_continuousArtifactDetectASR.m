function [winrej] = uf_continuousArtifactDetectASR(EEG, varargin)
%% bla
%   2. Same as 1. but per channel, using <locthresh>.
%
%
%
%   Input:
%     EEG:          continuous EEG dataset (EEGLAB's EEG structure)
%
%   Optional keywords:
%     'channels':   channels to be used for artifact detection. Can be a
%                   logical vector of length equal to EEG.chanlocs.
%                   Alternatively, numeric indeces, a cell with labels or
%                   a single regular expression.
%      'cutoff : Standard deviation cutoff for removal of bursts (via ASR). Data portions whose variance
%            is larger than this threshold relative to the calibration data are considered missing
%            data and will be removed. Default: 20 (following Chang 2019)
%   Output:
%     winrej:       winrej matrix flagging artifactual segments of data.
%                   Use with UF_CONITUOUSARTIFACTREJECT
%

cfg = finputcheck(varargin,...
    {'channel','integer',[],[];...
     'cutoff','real',[],20;...
    },'mode','error');
if ischar(cfg)
    error(cfg)
end
if isempty(cfg.channel)
    cfg.channel = 1:size(EEG.data,1);
end
% show an erro if eye-tracking channels are included
if isfield(EEG,'chanlocs') && ~isempty(EEG.chanlocs)
    assert(~any(~cellfun(@isempty, regexp({EEG.chanlocs(cfg.channel).labels},...
    'Eye|Pupil|EYE'))),...
    ['Detected Eye-Tracking channels. Usually including them is not a',...
    ' very good idea as they scale very different.']);
end

EEG.data = EEG.data(cfg.channel,:);
evalc("EEG = eeg_checkset(EEG);"); % evalc it to surpress warnings
evalc(sprintf("EEG_clean = clean_asr(EEG,%f,[],[],[],[],[],[],[],0);",cfg.cutoff)); % last argument activates riemann
% Following code adapted from "clean_artifcats.m" from 
% https://github.com/sccn/clean_rawdata
% detect what changed
sample_mask = ~(sum(abs(EEG.data-EEG_clean.data),1) < 1e-10);
% build winrej
winrej = reshape(find(diff([false sample_mask false])),2,[])';
winrej(:,2) = winrej(:,2)-1;

dur = sum(winrej(:, 2) - winrej(:, 1));

fprintf(['-----------------------------------------\n',...
    'ASR: Marked %.2f seconds as artifactual (%.2f%% of data)\n'...
    '-----------------------------------------\n'], dur / EEG.srate,...
    dur / EEG.pnts * 100);
end