function [resAll] = dc_spl2continuous(EEG,varargin)
%converts spline betas to continuous based values. This removes time-basis
% functions and parameter-basis functions
%
%Arguments:
%    'spline_idx': Which spline to evaluate?
%
%    'channel': Evaluate for specific channel (default all)
%
%    'time': Evaluate for specific points in time only (default all)
%
%    'spline_value': Specify at which values the spline should be evaluated
%
%    'deconv': Specifies whether to use dcBeta (default, deconv=1) or Xbeta
%      (the non-deconvoluted)
%
%Return:
%   beta (channels x times x 'predictors')
%
%*Example:*
% dc_spl2continuous(EEG,'spline_idx',1)



if isfield(EEG.deconv,'XBeta')
    nchan = size(EEG.deconv.XBeta,1);
elseif isfield(EEG.deconv,'dcBeta')
    nchan = size(EEG.deconv.dcBeta,1);
end


cfg = finputcheck(varargin,...
    {'spline_idx','integer',[],[]; %this is mandatory, which one do you want to evaluate?
    'channel','integer',[],1:nchan; % subselect channels
    'time','integer',[],[min(EEG.deconv.dcBasistime);max(EEG.deconv.dcBasistime)]; %take all time-points by default
    'spline_value','real',[],[]; %take all values defined in spline2valby default
    'deconv','boolean',[],1;
    },'mode','ignore');
if(ischar(cfg)); error(cfg);end


assert(cfg.spline_idx<=length(EEG.deconv.predictorSplines),'splineIdx could not be found in EEG.deconv')
assert(cfg.time(1)>=min(EEG.deconv.dcBasistime),'specified time is to small')
assert(cfg.time(2)<=max(EEG.deconv.dcBasistime),'specified time is to large')

cfg.time_selection = EEG.deconv.dcBasistime>=cfg.time(1) & EEG.deconv.dcBasistime<=cfg.time(2);
predictorSplines = EEG.deconv.predictorSplines{cfg.spline_idx};

if isempty(cfg.spline_value)
    cfg.spline_value_selection = 1:length(predictorSplines.spline2val);
else
    cfg.spline_value_selection = get_min(cfg.spline_value,predictorSplines.spline2val);
end

splineName = EEG.deconv.predictorSplines{cfg.spline_idx}.name;

betaIX = EEG.deconv.cols2variableNames ==  find(strcmp(EEG.deconv.variableNames,splineName));
% betaIX = ismember(EEG.deconv.colnames,predictorSplines.name);

if cfg.deconv
    b = EEG.deconv.dcBeta(cfg.channel,:,betaIX);
else
    b = EEG.deconv.XBeta(cfg.channel,:,betaIX);
end


% Could be optimized using mtimesx or maybe even gpuarray (but the overhead
% might be big?). It seems very fast anyway


% Push it in the samples-domain (get rid of time-splines)
tic

dat = nan(size(b,1),length(cfg.time_selection),size(b,3));
for chan = 1:size(b,1)
    if cfg.deconv
        dat(chan,:,:) =      EEG.deconv.dcBasis(:,cfg.time_selection)' * permute(b(chan,:,:),[2 3 1]); %permute instead of squeeze to keep order if chan = 1
    else
        dat(chan,:,:) = pinv(EEG.deconv.dcBasis(:,cfg.time_selection)) * permute(b(chan,:,:),[2 3 1]); %permute instead of squeeze to keep order if chan = 1
    end

end


% Push it in the value-domain (get rid of spline-predictors)
resAll = nan(size(b,1),size(EEG.deconv.dcBasis,2),length(cfg.spline_value_selection));
for chan = 1:size(b,1)
    resAll(chan,:,:) = permute(dat(chan,:,:),[2 3 1]) * predictorSplines.basis(cfg.spline_value_selection,:)'; %permute instead of squeeze to keep order if chan = 1
end
toc
