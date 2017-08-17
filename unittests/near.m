% Disclaimer: This function is taken from EEGlab and was NOT modified

% Testcases for EEGLab
% Copyright (C) 2006 Andreas Romeyke & Ronny Lindner
% Max-Planck-Institute for Human Cognitive and Brain Sciences Leipzig, Germany
% romeyke@cbs.mpg.de, art1@it-netservice.de
%
% based on EEGLab-toolbox
% http://www.sccn.ucsd.edu/eeglab/
% Copyright (C) 1996 Scott Makeig et al, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% this function checks if 2 arrays are equal using the epsilon value


function result = near(arg1, arg2, recursionlevel)
%      result = 0;
%  
%      if xor(isempty(arg1), isempty(arg2)) % 1 empty, the other not
%          return;
%      end;
%  
%      if isnumeric(arg1) && isnumeric(arg2)
%          if all(all( ...
%                        ((arg1-tc_epsilon <= arg2) & (arg2 <= arg1+tc_epsilon)) ...
%                      | (isnan(arg1) & isnan(arg2)) ...
%                    )) % workaround because Inf-Inf => NaN
%              result = 1;
%          end;
%      else
%          if isequal(arg1, arg2)
%              result = 1;
%          end;
%      end;

    result = 0;

    if nargin < 2
        return;
    end;

    if nargin < 3 || isempty(recursionlevel)
        recursionlevel = 5; % for structures
    end;

    if recursionlevel < 1
        return;
    end;

    if xor(isempty(arg1), isempty(arg2)) % 1 empty, the other not
        return;
    end;

    if isnumeric(arg1) && isnumeric(arg2)
        if ~isequal(size(arg1), size(arg2))
            return;
        end;
        if all(all( ...
                      ((arg1-tc_epsilon <= arg2) & (arg2 <= arg1+tc_epsilon)) ...
                    | (isnan(arg1) & isnan(arg2)) ...
                  )) % workaround because Inf-Inf => NaN
            result = 1;
        end;
    else
        if isstruct(arg1) && isstruct(arg2)
            f1 = fieldnames(arg1);
            f2 = fieldnames(arg2);
            if length(f1) == length(f2) && all(ismember(f1, f2))
                for num = 1:length(f1)
                    result = near([arg1.(f1{num})], [arg2.(f1{num})], recursionlevel-1);
                    if result == 0
                        return;
                    end;
                end;
            end;
        else
            if isequal(arg1, arg2)
                result = 1;
            end;
        end;
    end;

    % Testcases for EEGLab
% Copyright (C) 2006 Andreas Romeyke & Ronny Lindner
% Max-Planck-Institute for Human Cognitive and Brain Sciences Leipzig, Germany
% romeyke@cbs.mpg.de, art1@it-netservice.de
%
% based on EEGLab-toolbox
% http://www.sccn.ucsd.edu/eeglab/
% Copyright (C) 1996 Scott Makeig et al, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% this function returns the epsilon value
% used by floating point comparison in testcases
function tc_epsilon = tc_epsilon    
    tc_epsilon = 0.0001;
    