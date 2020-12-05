function [combi_winrej] = uf_combineWinrej(varargin)
%UF_COMBINE_WINREJ combines winrej arrays
%
%   This function can be used to combine the winrej arrays of multiple
%   artifact detection algorithms, before plugging the rejection array into
%   uf_continuousArtifactExclude.m
%
%   Input:
%           Any number of n*2 matrices with rejection latencies, as
%              produced by uf_continuousArtifactDetect() and
%              uf_continuousJointProbArtifactDetect().
%           Optionally, an EEG structure can be passed to this function,
%               to print information about the amount of data selected.
%
%   Example:
%           winrej_amp = uf_continuousArtifactDetect(EEG);
%           winrej_jp = uf_continuousJointProbArtifactDetect(EEG);
%           comb_winrej = uf_combine_winrej(winrej_amp, winrej_jp);
%           comb_winrej = uf_combine_winrej(winrej_amp, winrej_jp, EEG);
%           EEG = uf_continuousArtifactExclude(EEG, 'winrej', comb_winrej);
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

fprintf('\nCombining information from multiple artifact detection algorithms...');

% parse input
EEG_candidates = cellfun(@isstruct, varargin);
winrej_candidates = cellfun(@(x) isnumeric(x) & size(x, 2) == 2, varargin);
winrej_empty_candidates = cellfun(@(x) isnumeric(x) & isempty(x), varargin); % in case one of the algorithms did not detect any artifacts

% check assumptions
assert(all(EEG_candidates | winrej_candidates | winrej_empty_candidates),...
    'I don''t understand input nr. %i',...
    find(~(EEG_candidates | winrej_candidates | winrej_empty_candidates))) 
assert(sum(EEG_candidates) < 2, 'More than one EEG structure found');
winrejcheck = @(y) all(all(arrayfun(@(x) mod(x,1) == 0, y)));
assert(all(cellfun(winrejcheck, varargin(winrej_candidates))),...
    'winrej arrays should only contain latencies in integer indeces');

% if only one method detected artefacts, the function combines adjacent
% segments. The excluded samples are the same.

% print a message telling how many winrej arrays are combined
fprintf('\nCombining information from %i algorithms. \n', sum(winrej_candidates));

% combine info and sort by start latency
combi_winrej = cat(1, varargin{winrej_candidates});
combi_winrej = sortrows(combi_winrej);

% now find overlapping segments of artifactual data and combine them.
irow = 0;
while irow < size(combi_winrej, 1) - 1
    irow = irow + 1;
    this_end = combi_winrej(irow, 2);
    next_start = combi_winrej(irow + 1, 1);
    if next_start <= this_end
        combi_winrej(irow, 2) = combi_winrej(irow + 1, 2);
        combi_winrej(irow + 1, :) = [];
        % restart
        irow = 0;
    end
end

% print the optional summary stats:
if any(EEG_candidates)
    nbad = sum(combi_winrej(:, 2) - combi_winrej(:, 1));
    fprintf('\n %i segments were marked.', size(combi_winrej, 1));
    fprintf(['\n A total of %i samples (%.02f percent of data) was ',...
        'marked as bad.\n\n'], nbad,...
        nbad / size(varargin{EEG_candidates}.data, 2) * 100);
end

end

