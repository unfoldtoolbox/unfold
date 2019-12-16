function [Results,Info]= uf_tfce(cfg,D)
% data should be difference values with  Channels x Times x Subject x
% Condition


% Behinger version
% taken and adapted from ept_TFCE

% Licence from original:
% Copyright(C) 2012  Armand Mensen (14.12.2010)

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


tic
cfg = finputcheck(cfg,...
    {'E_H',   '', [], [0.66 2];...
    'nperm','integer',[],1000;...
    'neighbours','integer',[],[];...
    'chanlocs','struct',[],struct();...
    'statFun','',[],@limo_yuend_ttest;
    },'mode','ignore');


if(ischar(cfg)); error(cfg);end


%% Set Defaults

if ndims(D) == 2
    warning('found only two dimensions, I assume you have only a single channel')
    D = permute(D,[4 1 2 3]); %add an empty third dimension
end
if size(D,2) == 1
    disp('only a single channel found, setting neighbours to empty')
    cfg.neighbours = [];
else
    if isempty(cfg.neighbours)
        if isempty(cfg.chanlocs)
            error('please specify the chanlocs in your config')
        end
        
        disp('Calculating Channel Neighbours...')
        cfg.neighbours= ept_ChN2(cfg.chanlocs);
        disp('Done')
    end
end

%% Create all variables in loop at their maximum size to increase performance

maxTFCE = zeros(cfg.nperm,1);

%% Calculate the actual T values of all data

disp('Calculating Actual Differences...')

T_Obs = cfg.statFun(D(:,:,1:end-1,1),D(:,:,1:end-1,2));



% check for non-zero T-values (will crash TFCE)
if max(abs(T_Obs(:))) < 0.00001
    error('T-values were all 0')
end

TFCE_Obs = ept_mex_TFCE2D(T_Obs, cfg.neighbours, cfg.E_H);

%% Calculating the T value and TFCE enhancement of each different permutation

disp('Calculating Permutations...')

for j   = 1:cfg.nperm           % Look into parfor for parallel computing
    
    
    Signs =[-1,1];
    SignSwitch = randsample(Signs,size(D,3),'true')';
    
    Dperm = permute(SignSwitch.*permute(D,[3 1 2 4]),[2 3 1 4]);
    [T_Perm,~] = cfg.statFun(Dperm(:,:,:,1),Dperm(:,:,:,2));
    
    
    TFCE_Perm = ept_mex_TFCE2D(T_Perm, cfg.neighbours, cfg.E_H);
        
    maxTFCE(j) = max(abs(TFCE_Perm(:)));       % stores the maximum absolute value
    
    progressbar(j/cfg.nperm); %Progress Bar
    
end

disp('Done')

%% Calculating the p value from the permutation distribution

disp('Calculating P-Values and Saving...')

% add observed maximum
edges = [maxTFCE;max(abs(TFCE_Obs(:)))];

[~,bin]     = histc(abs(TFCE_Obs),sort(edges));
P_Values    = 1-bin./(cfg.nperm+2);

% Save test information in single structure
c = clock;
Info.Comments = ['TFCE analysis conducted at ', num2str(c(4)), ':', num2str(c(5)), ' on ', date];

Info.Parameters.E_H         = cfg.E_H;
Info.Parameters.nperm       = cfg.nperm;

Info.Electrodes.ChannelNeighbours = cfg.neighbours;

Results.Obs                 = T_Obs;
Results.TFCE_Obs            = TFCE_Obs;
Results.maxTFCE             = sort(maxTFCE);
Results.P_Values            = P_Values;

%%
disp('All done!')
toc

[min_P, idx] = min(Results.P_Values(:));
[Ch, S]      = ind2sub(size(Results.P_Values),idx);
max_Obs      = Results.Obs(idx);

display(['Peak significance found at channel ', num2str(Ch), ' at sample ', num2str(S), ': max-t-val(', num2str(size(D,1)-1), ') = ', num2str(max_Obs), ', p = ', num2str(min_P)]);


end
