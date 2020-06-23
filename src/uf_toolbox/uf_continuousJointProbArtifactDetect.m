function [winrej] = uf_continuousJointProbArtifactDetect(EEG, varargin)
%% winrej = UF_CONTINUOUSJOINTPROBARTIFACTDETECT(EEG, varargin)
%
%   Mark improbable data segments in continuous data. Intended for use with
%   the unfold-toolbox (github.com/unfoldtoolbox/unfold). Can be used
%   instead of UF_CONTINUOUSARTIFACTDETECT or in addition to it.
%   This function is based on eeglab's JOINTPROB function.
%
%   Improbable data are marked in two ways:
%   1. A probability distribution using all channels across segments is
%      created and a segment is marked as "improbable", when all channels
%      deviate by <globthresh> * std.dev.
%      As this function works on continuous data, pseudo-segments of length
%      <seglength> are created.
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
%     'locthresh':  activity probability limit in std.dev.; either a single
%                   value applied to all channels or a vector with one
%                   value per channel.
%     'globthresh': global limit(s) (all activities grouped) (in std. dev.)
%     'seglength':  data need to be pseudo-segmented for joint-probability.
%                   How long these pseudo-segments are also defines how
%                   long the marked artifactual periods are. Provide a
%                   value in seconds (default 0.5).
%     'verbose':    Display messages created by 'jointprob'? (default: false)
%
%   Output:
%     winrej:       winrej matrix flagging artifactual segments of data.
%                   Use with UF_CONITUOUSARTIFACTREJECT
%
%   Author: Wanja Moessing, w.a.moessing@gmail.com, 13/05/2020
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


cfg = finputcheck(varargin,...
    {'globthresh',   'real', [], 3;...
    'locthresh',     'real', [], 3;...
    'robust_normalization','boolean',[],1;...
    'seglength','real',[],0.5;...
    'verbose','boolean',[],0;...
    'channel','integer',[],[];...
    'plot','boolean',[],0;
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

%% joint probability is calculated over segments of data, so we need to
% create artificial epochs including the whole continuous data.
% Removing any baseline other than the average does not make sense, as we
% wouldn't know what's happening there. Not removing the baseline might
% harm joint probability.

% some rare imprecisions in EEG.xmax lead to incomplete segmentation
rounderr = (1 - EEG.xmax / (EEG.pnts / EEG.srate));

% check if such an imprecision is present and if it's just one sampling
% point. If so, simply correct that.
if EEG.xmax ~= EEG.pnts / EEG.srate &&...
        round(rounderr * EEG.pnts) <= 1
    fprintf('Correcting round-off error of %.6f seconds in EEG.xmax\n',...
        rounderr);
    EEG.xmax = EEG.pnts/EEG.srate;
end

% remove all event fields that aren't necessary for joint probability.
% eeg_checkset in eeg_regepoch will otherwise take ages to make them
% uniform and the temporary dataset will not be used after this script,
% anyways.
fnames = fieldnames(EEG.event);
fnames = fnames(~ismember(fnames, {'type', 'latency', 'epoch'}));
EEG.event = rmfield(EEG.event, fnames);

disp('creating pseudo-epochs...');
% this also creates pseudo events, which we later use to get the latencies
if cfg.verbose
    EEGreg = eeg_regepochs(EEG, 'recurrence', cfg.seglength,...
        'rmbase', -cfg.seglength, 'limits', [0, cfg.seglength],...
        'eventtype', 'JP_artificial');
else
    evalc(['EEGreg = eeg_regepochs(EEG, ''recurrence'', cfg.seglength,',...
        '''rmbase'', -cfg.seglength, ''limits'', [0, cfg.seglength],',...
        '''eventtype'', ''JP_artificial'');']);
end

% check data
assert(((EEGreg.pnts * EEGreg.trials) / EEG.pnts) == 1,...
    'JP: segmented length of data does not equal continuous length!');


%% run jointprob
%EEG = pop_jointprob(EEG, 1, 1:EEG.nbchan, P.locthresh,...
%    P.globthresh, 1, 0, 1, [], P.plot);

% taken from pop_jointprob

% local
[~,rejEtmp ] = jointprob( EEGreg.data(cfg.channel,:,:), cfg.locthresh,[], cfg.robust_normalization+1); % +1 because 0 = no-norm, 1 = norm, 2=trimmed norm

% "global", i.e. concatenating over channels
tmpdata2 = permute(EEGreg.data(cfg.channel,:,:), [3 1 2]);
tmpdata2 = reshape(tmpdata2, size(tmpdata2,1), size(tmpdata2,2)*size(tmpdata2,3));
[tmp,rejG ] = jointprob( tmpdata2, cfg.globthresh, [], cfg.robust_normalization+1); 

rejE    = zeros(size(EEGreg.data,1), size(rejEtmp,2));
rejE(cfg.channel,:) = rejEtmp;
rej = rejG' | max(rejE, [], 1);


if cfg.plot
rejstatepoch( EEGreg.data(cfg.channel,:,:), rejEtmp, 'global', 'on', 'rejglob', rejG, ...
						'threshold', cfg.locthresh, 'thresholdg', cfg.globthresh, 'normalize', 'on')
end
%% translate to winrej for uf_continuousArtifactExclude
% for artifactual trials, identify the indeces of the artificial triggers created above
badtrls = {EEGreg.epoch(rej).eventtype};

% by definition, we want the first and last match per "epoch"
uridx = cellfun(@(x) quantile(find(strcmp(x, 'JP_artificial')), [0, 1]),...
    badtrls, 'UniformOutput', 0);

% extract all urevents
urevents = {EEGreg.epoch(rej).eventurevent};

% extract latencies
winrej = zeros(length(uridx), 2);
for ilat = 1: length(uridx)
    foo = [urevents{ilat}{uridx{ilat}}];
    winrej(ilat, 1:2) = [EEGreg.urevent(foo).latency];
end

% display results: (% should match what uf_continuousArtifactExclude says)
dur = sum(winrej(:, 2) - winrej(:, 1));
durE = sum(max(rejE, [], 1))*cfg.seglength*EEG.srate;
durG = sum(rejG)*cfg.seglength*EEG.srate;

fprintf(['-----------------------------------------\n',...
    'JointProbability: Marked %.2f seconds as artifactual (%.2f%% of data)\n'...
    '%.2f seconds due to local alone, %.2f seconds due to global alone \n',...
    '-----------------------------------------\n'], dur / EEGreg.srate,...
    dur / EEG.pnts * 100,durE/EEGreg.srate,durG/EEGreg.srate);
end