function ufresult = uf_addmarginal(ufresult,varargin)
%add the marginal of the other predictors (i.e. continuous & spline
%predictors) to the beta estimates.
% Important: If dummy-coded (i.e. non-effect coded) predictors and
% interactions exist, they are NOT added to the marginal effect. I.e. the
% output of the method returns the average ERP evaluated at the average of
% all spline/continuous predictors, keeping the categorical/interaction
% structure untouched.
%
%
% For instance the model 1 + cat(facA) + continuousB
% has the betas: intercept, facA==1, continuousB-Slope
%
% The beta output of uf_condense(uf_glmfit) mean the following:
% intercept: response with facA = 0 and continuousB = 0
% facA==1  : differential effect of facA == 1 (against facA==0)
% continuousB-slope: the slope of continous B
%
% Using uf_predictContinuous, we evaluate the continuous parameter at [0 50 100]
% The beta output of uf_predictContinuous mean the following:
% intercept: same as before
% facA==1  : same as before
% continuousB@0  : the differential effect if continuous B is 0 
% continuousB@50 : the differential effect if continuous B is 50
% continuousB@100: the differential effect if continuous B is 100
%
% Using uf_addmarginal, the average response is added to all predictors.
%
% intercept: the response of facA==0 AND continuousB@mean(continuousB)
% intercept: the response of facA==1 AND continuousB@mean(continuousB)
% continuousB@0  : the response of facA==0 if continuous B is 0 
% continuousB@50 : the response of facA==0 if continuous B is 50
% continuousB@100: the response of facA==0 if continuous B is 100
%
% Note that mean(continuousB) does not need to be a number we evaluated in
% the uf_predictContinuous step.



% parse inputs
cfg = finputcheck(varargin,...
    {'channel','integer',[],[]; ...
    'betaSetname','string','','' ...
    },'mode','ignore');

if(ischar(cfg)); error(cfg);end

% In order to add the marginal, we need evaluated splines (uf_predictContinuous) first.
% here we are looking for spline_converted or continuous_converted. I.e. if
% any have not been converted, throw an error.
if any(strcmp({ufresult.param.type},'spline')) || any(strcmp({ufresult.param.type},'continous'))
    error('In order to add the marginals, you need to run uf_predictContinuous first to evaluate the splines and continuous predictors');
end


if isempty(cfg.betaSetname)
    [betaSetname] = uf_unfoldbetaSetname(ufresult,varargin{:});
    
    % RECURSION ALERT!
    if length(betaSetname)>1
        for b = betaSetname
            ufresult_tmp = uf_addmarginal(ufresult,'betaSetname',b{1});
            ufresult.(b{1}) = ufresult_tmp.(b{1});
        end
        return
    else
        cfg.betaSetname = betaSetname{1};
    end
end
% END OF RECURSION ALERT



fprintf('uf_addmarginal: working on field %s \n',cfg.betaSetname)
if isempty(cfg.channel)
    cfg.channel = 1:size(ufresult.(cfg.betaSetname),1);
end
event = [ufresult.param.event];
paramList = {ufresult.param.name};

fprintf('re-running uf_condense to recover unconverted splines \n')
ufresult_avg = uf_condense(ufresult);
ufresult_avg = uf_predictContinuous(ufresult_avg,'auto_method','average');
for e = unique(event)
    % check if predictor and event combination exists
    eventIdx = find(strcmp(e,event));
    if length(eventIdx) == 1
        continue
    end
    eventParam = paramList(eventIdx);
    for p = unique(eventParam)
        % Find the names & types of the other parameters
        currEvent = eventIdx(strcmp(p,eventParam));
        otherEvents = setdiff(eventIdx,currEvent);
        otherParamNames = unique(paramList(otherEvents));
        ufresultavg_ix = [];
        for pOther = otherParamNames
            ufresultavg_ix(end+1) = find(strcmp(pOther{1},{ufresult_avg.param.name}));
        end
        % exclude all categoricals.
        % they should be covered by the effects/dummy coding schema
        % already ?!
        removeix = strcmp('categorical',{ufresult_avg.param(ufresultavg_ix).type});
        removeix = removeix|strcmp('interaction',{ufresult_avg.param(ufresultavg_ix).type});
        ufresultavg_ix(removeix) = [];
        % calculate the marginal over all other predictors
        average_otherEffects = squeeze(sum(ufresult_avg.(cfg.betaSetname)(cfg.channel,:,ufresultavg_ix),3));
        
        % add this marginal to the current predictor
        ufresult.(cfg.betaSetname)(cfg.channel,:,currEvent) = ufresult.(cfg.betaSetname)(cfg.channel,:,currEvent) + repmat(average_otherEffects,1,1,length(currEvent));
        
    end
end