function unfold = dc_addmarginal(unfold,varargin)
%add the marginal of the other predictors (i.e. continuous & spline
%predictors to the beta estimates

% parse inputs
cfg = finputcheck(varargin,...
    {'channel','integer',[],[]; ...
    'betaSetname','string','','' ...
    },'mode','ignore');

if(ischar(cfg)); error(cfg);end

% In order to add the marginal, we need evaluated splines (dc_getParam) first.
% here we are looking for spline_converted or continuous_converted. I.e. if
% any have not been converted, throw an error.
if any(strcmp({unfold.param.type},'spline')) || any(strcmp({unfold.param.type},'continous'))
    error('In order to add the marginals, you need to run dc_getParam first to evaluate the splines and continuous predictors');
end


if isempty(cfg.betaSetname)
    [betaSetname] = dc_unfoldbetaSetname(unfold,varargin{:});
    
    % RECURSION ALERT!
    if length(betaSetname)>1
        for b = betaSetname
            unfold_tmp = dc_addmarginal(unfold,'betaSetname',b{1});
            unfold.(b{1}) = unfold_tmp.(b{1});
        end
        return
    else
        cfg.betaSetname = betaSetname{1};
    end
end
% END OF RECURSION ALERT



fprintf('dc_addmarginal: working on field %s \n',cfg.betaSetname)
if isempty(cfg.channel)
    cfg.channel = 1:size(unfold.(cfg.betaSetname),1);
end
event = [unfold.param.event];
paramList = {unfold.param.name};

fprintf('re-running beta2unfold to recover unconverted splines \n')
unfold_avg = dc_beta2unfold(unfold);
unfold_avg = dc_getParam(unfold_avg,'auto_method','average');
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
        unfoldavg_ix = [];
        for pOther = otherParamNames
            unfoldavg_ix(end+1) = find(strcmp(pOther{1},{unfold_avg.param.name}));
        end
        % exclude all categoricals.
        % they should be covered by the effects/dummy coding schema
        % already ?!
        removeix = strcmp('categorical',{unfold_avg.param(unfoldavg_ix).type});
        removeix = [removeix strcmp('interaction',{unfold_avg.param(unfoldavg_ix).type})];
        unfoldavg_ix(removeix) = [];
        % calculate the marginal over all other predictors
        average_otherEffects = squeeze(sum(unfold_avg.(cfg.betaSetname)(cfg.channel,:,unfoldavg_ix),3));
        
        % add this marginal to the current predictor
        unfold.(cfg.betaSetname)(cfg.channel,:,currEvent) = unfold.(cfg.betaSetname)(cfg.channel,:,currEvent) + repmat(average_otherEffects,1,1,length(currEvent));
        
    end
end