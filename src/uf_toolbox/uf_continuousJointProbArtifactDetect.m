function [winrej] = uf_continuousJointProbArtifactDetect(EEG, globthresh, locthresh, varargin)
%winrej = UF_CONTINUOUSJOINTPROBARTIFACTDETECT(EEG, varargin)
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

% parse variable input
p = inputParser;
p.FunctionName = 'uf_continuousJointProbArtifactDetect';
p.addRequired('EEG', @isstruct);
p.addRequired('globthresh', @(x) isnumeric(x) && isscalar(x));
p.addRequired('locthresh', @(x) isnumeric(x) && isscalar(x));
p.addParameter('seglength', 0.5, @(x) isnumeric(x) && isscalar(x));
p.addParameter('verbose', false, @islogical);
p.addParameter('channels',':',...
    @(x) (ischar(x) | islogical(x) | isnumeric(x) | iscell(x)));
p.addParameter('plot', 0, @(x) islogical(x) || ismembc(x, [0,1]));
parse(p, EEG, globthresh, locthresh, varargin{:})
P = p.Results;

%% parse channel input
if strcmp(P.channels, ':')
    P.channels = 1:length({EEG.chanlocs.labels});
else
    if all(size(P.channels) == 1)
        if iscell(P.channels)
            P.channels = cell2mat(P.channels);
        end
        switch class(P.channels)
            case 'char' % regex/single channel label
                chantype = 'regex';
            case 'double' % single channel index
                chantype = 'index';
        end
    else % logical, numeric indeces, label vector
        % case cell vector
        if iscell(P.channels)
            %sanity check: all inputs the same class?
            classes = cellfun(@metaclass, P.channels, 'UniformOutput', 0);
            assert(all(strcmp(classes, classes{1})),...
                'Please specify channels EITHER as numbers OR as labels!');
            if strcmp(classes{1}, 'char')
                chantype = 'labels';
            elseif strcmp(classes{1}, 'logical')
                P.channels = cell2mat(P.channels);
                chantype = 'logindex';
            elseif isnumeric(P.channels{1})
                P.channels = cell2mat(P.channels);
                if all(ismembc(P.channels, [0, 1]))
                    chantype = 'logindex';
                else
                    chantype = 'index';
                end
            end
        elseif isnumeric(P.channels) % case numeric vector
            if all(ismembc(P.channels, [0, 1]))
                chantype = 'logindex';
            else
                chantype = 'index';
            end
        elseif islogical(P.channels) % case logic indeces
            chantype = 'logindex';
        elseif ischar(P.channels)
            chantype = 'regex';
        end
    end
    
    % pop_select can handle labels and indexes
    % in case of regex, first convert that to a logical index
    if strcmp(chantype, 'regex')
        P.channels = ~cellfun(@isempty, regexp({EEG.chanlocs.labels},...
            P.channels));
    end
    
    % convert logical to indeces
    if ismember(chantype, {'regex', 'logindex'})
        P.channels = find(P.channels);
    end
end

% show an erro if eye-tracking channels are included
assert(~any(~cellfun(@isempty, regexp({EEG.chanlocs(P.channels).labels},...
    'Eye|Pupil|EYE'))),...
    ['Detected Eye-Tracking channels. Usually including them is not a',...
    ' very good idea as they scale very different.']);

%% create a reduced set without eyetracking etc.
fprintf(['\nuf_continuousJointProbArtifactDetect: the following changes '...
    'are only done internally to detect artefactual latency ranges..\n'])
JP_EEG = pop_select(EEG, 'channel', P.channels);

%% joint probability is calculated over segments of data, so we need to
% create artificial epochs including the whole continuous data.
% Removing any baseline other than the average does not make sense, as we
% wouldn't know what's happening there. Not removing the baseline might
% harm joint probability.

% some rare imprecisions in EEG.xmax lead to incomplete segmentation
rounderr = (1 - JP_EEG.xmax / (JP_EEG.pnts / JP_EEG.srate));

% check if such an imprecision is present and if it's just one sampling
% point. If so, simply correct that.
if JP_EEG.xmax ~= JP_EEG.pnts / JP_EEG.srate &&...
        round(rounderr * JP_EEG.pnts) <= 1
    fprintf('Correcting round-off error of %.6f seconds in EEG.xmax\n',...
        rounderr);
    JP_EEG.xmax = JP_EEG.pnts/JP_EEG.srate;
end

% remove all event fields that aren't necessary for joint probability.
% eeg_checkset in eeg_regepoch will otherwise take ages to make them
% uniform and the temporary dataset will not be used after this script,
% anyways.
fnames = fieldnames(JP_EEG.event);
fnames = fnames(~ismember(fnames, {'type', 'latency', 'epoch'}));
JP_EEG.event = rmfield(JP_EEG.event, fnames);

disp('creating pseudo-epochs...');
% this also creates pseudo events, which we later use to get the latencies
if P.verbose
    JP_EEG = eeg_regepochs(JP_EEG, 'recurrence', P.seglength,...
        'rmbase', -P.seglength, 'limits', [0, P.seglength],...
        'eventtype', 'JP_artificial');
else
    evalc(['JP_EEG = eeg_regepochs(JP_EEG, ''recurrence'', P.seglength,',...
        '''rmbase'', -P.seglength, ''limits'', [0, P.seglength],',...
        '''eventtype'', ''JP_artificial'');']);
end
% check data
assert(((JP_EEG.pnts * JP_EEG.trials) / EEG.pnts) == 1,...
    'JP: segmented length of data does not equal continuous length!');

%% run jointprob
JP_EEG = pop_jointprob(JP_EEG, 1, 1:JP_EEG.nbchan, P.locthresh,...
    P.globthresh, 1, 0, 1, [], P.plot);

%% translate to winrej for uf_continuousArtifactExclude
% for artifactual trials, identify the indeces of the artificial triggers created above
badtrls = {JP_EEG.epoch(JP_EEG.reject.rejjp).eventtype};

% by definition, we want the first and last match per "epoch"
uridx = cellfun(@(x) quantile(find(strcmp(x, 'JP_artificial')), [0, 1]),...
    badtrls, 'UniformOutput', 0);

% extract all urevents
urevents = {JP_EEG.epoch(JP_EEG.reject.rejjp).eventurevent};

% extract latencies
winrej = zeros(length(uridx), 2);
for ilat = 1: length(uridx)
    foo = [urevents{ilat}{uridx{ilat}}];
    winrej(ilat, 1:2) = [JP_EEG.urevent(foo).latency];
end

% display results: (% should match what uf_continuousArtifactExclude says)
dur = sum(winrej(:, 2) - winrej(:, 1));
fprintf(['-----------------------------------------\n',...
    'JointProbability: Marked %.2f seconds as artifactual ',...
    '(%.2f%% of data)\n',...
    '-----------------------------------------\n'], dur / JP_EEG.pnts,...
    dur / EEG.pnts * 100);
end