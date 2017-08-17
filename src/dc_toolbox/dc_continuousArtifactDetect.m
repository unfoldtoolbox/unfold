function [WinRej] = dc_continuousArtifactDetect(EEG,varargin)
%Reject commonly recorded artifactual potentials (c.r.a.p.)
% This function has been altered very much by Benedikt Ehinger
% I removed all the filter-features
% I changed the input parser.
%
% There are a number of common artifacts that you will see in nearly every EEG data file. These
% include eyeblinks, slow voltage changes (caused mostly by skin potentials), muscle activity
% (from moving the head or tensing up the muscles in the face or neck), horizontal eye movements,
% and various types of C.R.A.P. (Commonly Recorded Artifactual Potentials).
%
% Although we usually perform artifact rejection on the segmented data, it's a good idea to
% examine the raw unsegmented EEG data first. You can usually identify patterns of artifacts,
% make sure there were no errors in the file, etc., more easily with the raw data [1].
%
% crap.m allows you to automatically identify large peak-to-peak differences or extreme amplitude
% values, within a moving window, across your continuous EEG dataset. After performing crap.m,
% artifactual segments will be rejected and replaced by a 'boundary' event code.
%
%Arguments:
% EEG:         - continuous EEG dataset (EEGLAB's EEG structure)
% 'amplitudeThreshold':     - Thresolds ( values). [-lim +lim] is marked
% 'windowsize':     - moving window width in msec (default 2000 ms)
% 'stepsize':    - moving window step (default 1000 ms)
% 'combineSegements':  - marked segment(s) closer than this value will be joined together.
%
%Reference:
% ERP Boot Camp: Data Analysis Tutorials. Emily S. Kappenman, Marissa L. Gamble, and Steven J. Luck. UC Davis
%
% This function is part of ERPLAB Toolbox
% Author: Javier Lopez-Calderon
% Center for Mind and Brain
% University of California, Davis,
% Davis, CA
% 2009

'';
%
% ERPLAB Toolbox
% Copyright � 2007 The Regents of the University of California
% Created by Javier Lopez-Calderon and Steven Luck
% Center for Mind and Brain, University of California, Davis,
% javlopez@ucdavis.edu, sjluck@ucdavis.edu
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Acknowledgment:
% Thanks to Paul Kieffaber for his advice on improving this function.
%


cfg = finputcheck(varargin,...
    {'amplitudeThreshold',   'integer', [], 150;...
    'windowsize','integer',[],2000;...
    'channels','integer',[1:size(EEG.data,1)],1:size(EEG.data,1);
    'stepsize','integer',[],100;...
    'combineSegements','integer',[],[];...
    },'mode','ignore');


if(ischar(cfg)); error(cfg);end


ampth     = cfg.amplitudeThreshold;%microV
winms     = cfg.windowsize;%ms
stepms    = cfg.stepsize;% ms steps
chanArray = cfg.channels;


shortisi  = cfg.combineSegements;



[WinRej, chanrej] = basicrap(EEG, chanArray, ampth, winms, stepms);
shortisisam  = floor(shortisi*EEG.srate/1000);  % to samples

if isempty(WinRej)
    fprintf('\n No data found worse than the threshold. No rejection was performed.\n');
else
    %colorseg = [1.0000    0.9765    0.5294];
    if ~isempty(shortisisam)
        [WinRej, chanrej ] = joinclosesegments(WinRej, chanrej, shortisisam);
    end
    % Correct overlapping windows, there used to be a bug in EEGLAB and this fixes it.
    throw_out = nan;
    while ~isempty(throw_out)
        throw_out = [];
        for i = 1:size(WinRej,1)-1
            if ~isempty(throw_out) && throw_out(end) == i
                continue
            end
            if WinRej(i,2)>=WinRej(i+1,1)
                throw_out=[throw_out i+1];
                WinRej(i,1) = min(WinRej(i,1),WinRej(i+1,1));
                WinRej(i,2) = max(WinRej(i,2),WinRej(i+1,2));
            end



        end
        WinRej(throw_out,:) = [];
    end


    % colormatrej = repmat([1 0 0], size(WinRej,1),1);
    % matrixrej = [WinRej colormatrej chanrej];
    % % call figure
    % eegplot(EEG.data, 'winrej', matrixrej, 'srate', EEG.srate,'events', EEG.event,'winlength', 50,'dispchans',1,'command','');
    fprintf('\n %g segments were marked.\n\n', size(WinRej,1));


end
