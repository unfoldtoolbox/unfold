function output = dc_getParam(unfold,varargin)
%% Evaluates a continuous/spline parameter at specific values
% For example you get an beta for par1 at 3 for a continuous
% variable. The output then would be the respective values 30,60 and
% 90. Because model-estimates / parameters are defined for each time-point&electrode and can
% also encompass multiple betas (in the case of splines), this
% becomes non trivial and thus this function.
%
%Arguments:
%   cfg.pred_value(cell): One entry per parameter:
%       {{'par1',[10 20 30]},{'par2',[0,1,2]}}.
%       This evaluates parameter 1 at the values 10,20 and 30. Parameter 2
%       at 0, 1 and 2.
%       Default behaviour: evaluates 7 linearly spaced values between the min + max. of the
%       parameterdomain
%   cfg.auto_method(string): 'quantile' (default) or 'linear'.
%       'quantile' - the auto_n values are placed on the quantile of the predictor
%       'linear'   - the auto_n
%   cfg.auto_n (integer) : default 7; the number of automatically evaluated values 
%
%Return:
%   Result-Betas with evaluated betas at specified continuous values.
%
%*Example:*
% TODO


cfg = finputcheck(varargin,...
    {'pred_value','cell',[],{{'',[]}};
    'deconv','integer',[-1,0,1],-1;
    'convertSplines','integer',[],1;
    'auto_method','string',{'quantile','linear'},'quantile';
    'auto_n','integer',[],10;
    },'mode','ignore');

if(ischar(cfg)); error(cfg);end


beta_dcExists   = isfield(unfold,'beta')&& isnumeric(unfold.beta);
beta_nodcExists = isfield(unfold,'beta_nodc') && isnumeric(unfold.beta_nodc);

if cfg.deconv == 1
    assert(beta_dcExists,'beta_dc missing or not numeric')
    
    
elseif cfg.deconv == 0
    assert(beta_nodcExists,'beta_nodc missing or not numeric')
    
    
elseif cfg.deconv == -1 % auto detect, recursive call
    
    assert(isfield(unfold,'beta')|isfield(unfold,'beta_nodc'),'error: to use autodetect at least the field unfold.beta  or unfold.beta_nodc needs to exist')
    fn = fieldnames(unfold);
    
    if isfield(unfold,'beta')
        sizeBeta = size(unfold.beta);
    else
        sizeBeta = size(unfold.beta_nodc);
    end
    %-------------- Recursive part
    
    for f = fn'
        if strcmp(f,'times')
            continue
        end
        if length(sizeBeta) == length(size(unfold.(f{1}))) &&  all(sizeBeta == size(unfold.(f{1})))
            
            % overwrite the "beta" field and use that
            cfg.deconv = 1;
            unfold_tmp = unfold;
            unfold_tmp.beta = unfold.(f{1});
            output_tmp = dc_getParam(unfold_tmp,cfg); %------- recursive call
            if ~exist('output','var')
                output = output_tmp;
                output = rmfield(output,'beta'); %delete the temporary beta field
            end
            output.(f{1}) = output_tmp.beta;
        end
    end
    
    
    return
    %-------------- End Recursive part
end


% Array of the sorts: {{'parName',linspace(0,10,5)},{'parname2',1:5}}
predValueSelectList = cfg.pred_value;
predNameList = cellfun(@(x)x{1},predValueSelectList,'UniformOutput',0);
% if cfg.convertSplines == 1
%      unfold = dc_beta2unfold(unfold,'convertSplines',1);
% end
if cfg.convertSplines == 1
    [splIdxList,paramList] = dc_getSplineidx(unfold);
    paramListSpline = paramList;
    
    splineList = strcmp(unfold.deconv.variableType,'spline');
    for spl = [1:sum(splineList);find(splineList)]
        ix = find(spl(2)== unfold.deconv.cols2variableNames,1,'first');
        ix = find(ix == paramList);
        NsplineVal = length(unfold.deconv.predictorSplines{spl(1)}.spline2val);
        NsplinePred = size(unfold.deconv.predictorSplines{spl(1)}.X,2);
        paramListSpline(ix+1:end) = paramListSpline(ix+1:end)+ NsplineVal - NsplinePred;
%         paramListSpline
    end
    
else
    paramList = 1: size(unfold.deconv.X,2);
    paramListSpline = paramList;
end
if cfg.deconv == 1
    beta = unfold.beta;
elseif cfg.deconv == 0
    
    beta = unfold.beta_nodc;
end

betaNew = [];
epochNew = unfold.epoch(1); %needs to be removed
for predIDX = [paramList;paramListSpline]
    variableIdx = unfold.deconv.cols2variableNames(predIDX(1));
    
    b = beta(:,:,predIDX(2));
    
    e = unfold.epoch(predIDX(2));
    
    if cfg.convertSplines && variableIdx ~=0 && strcmp(unfold.deconv.variableType{variableIdx},'spline')
        
        splName = cellfun(@(x)x.name,unfold.deconv.predictorSplines,'UniformOutput',0);
        splIdx = find(strcmp(splName,unfold.deconv.variableNames{variableIdx}));
        
        spl = unfold.deconv.predictorSplines{splIdx};
        % If we found the spline_value given value in the splines, then choose
        % that one
        customSplineValue = strcmp(predNameList,spl.name);
        if any(customSplineValue)
            splValueSelect = predValueSelectList{customSplineValue}{2};
            
        else
            splValueSelect = auto_spacing(cfg,spl.paramValues);
        end
        
        
        % all betas of the spline
        ix = strcmp({unfold.epoch.name},spl.name);
        if sum(ix) ~= length(spl.spline2val)
            error('something went wrong, maybe you did not convert the splines (''convertSplines'') while calling beta2unfold, or you ran dc_getParam already?')
        end
        
        indexList = find(ix);
        for c = splValueSelect
            thisBetaIx = get_min(c,[unfold.epoch(ix).value]);
            thisBetaIx = indexList(thisBetaIx);
            
            betaNew(:,:,end+1) =beta(:,:,thisBetaIx);
            eNew = unfold.epoch(thisBetaIx);
%             eNew.value = c;
%             eNew.name = spl.name;
%             eNew.event = unfold.deconv.eventtype{variableIdx};
            epochNew(end+1) = eNew;
        end
        
        %         t(strcmp(t.parameter,spl.name)& ~ismember(t.paramvalue,,:) = [];
    elseif variableIdx ~=0 && strcmp(unfold.deconv.variableType{variableIdx},'continuous')
        % we have either a categorical or a continuous predictor here
        customContValue = strcmp(predNameList,unfold.deconv.colnames(predIDX(1)));
        if any(customContValue)
            contValueSelect = predValueSelectList{customContValue}{2};
            
        else
            values = unfold.deconv.X(:,predIDX(1));
            % we do not have access at this point to which 'X' entry belongs to that event. 
            % Thus we discard every 0 value and hope a bit, that it does not matter because 0 is the intercept anyway.
            % This is supoptimal and I'm sorry if it creats inconveniences.
            % We would need to introduce a whole new field to carry around
            % to compensate for this.
            contValueSelect = auto_spacing(cfg,values(values~=0));
        end
        
        for c = contValueSelect
            betaNew(:,:,end+1) = b.*c;
            eNew = e;
            eNew.value = c;
            epochNew(end+1) = eNew;
        end
        
    else % categorical & or not converting
        betaNew(:,:,end+1) = b;
        
        epochNew(end+1) = e;
    end
    
end
 epochNew(1) = [];
betaNew(:,:,1) = [];

if cfg.deconv
    unfold.beta = betaNew;
else
    unfold.beta_nodc = betaNew;
end
unfold.epoch = epochNew;
output = unfold;

end

function contValueSelect = auto_spacing(cfg,predVal)
switch cfg.auto_method
    case 'linear'
        % I experienced a small bug where when I took the min/max values
        % directly I found either very noisy, or flat values. So i thought
        % it would be good to move 5% more to the mean
        if isempty(predVal) || length(unique(predVal)) == 1
            warning('skipping auto-variableselection')
            contValueSelect = 1;
            return
        end
        contmax = max(predVal);
        contmin = min(predVal);
        ran = contmax-contmin;
        
        contValueSelect = linspace(contmin+0.05*ran,contmax-0.05*ran,cfg.auto_n);
    case 'quantile'
        contValueSelect = quantile(predVal,linspace(0,1,cfg.auto_n));
end
end
