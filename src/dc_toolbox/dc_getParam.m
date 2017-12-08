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
%       'linear'   - the auto_n values are placed linearly over the range of the predictor
%       'average'  - only evaluates at the average of the predictor
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
    'auto_method','string',{'quantile','linear','average'},'quantile';
    'auto_n','integer',[],10;
    },'mode','ignore');

if(ischar(cfg)); error(cfg);end

% check if function has been run before
% this will lead to mistakes and errors because some betas will not be the betas we expect anymore - don't allow that!
if any(cellfun(@(x)~isempty(x),strfind({unfold.param.type},'converted')))
    error('cannot run dc_getParam twice. Run dc_beta2unfold again first')
end


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
[~,paramList] = dc_getSplineidx(unfold);
if cfg.deconv == 1
    beta = unfold.beta;
elseif cfg.deconv == 0
    
    beta = unfold.beta_nodc;
end

betaNew = [];
epochNew = unfold.param(1); %needs to be removed
for currPred= 1:length(paramList)
    predIDX = paramList(currPred);
    
    % if we are at the last predictor, this one goes to the end of the
    % designmatrix
    if currPred == length(paramList)
        predIDX_next = size(unfold.deconv.X,2);
    else
        %else it goes to the next predictor
        predIDX_next = paramList(currPred+1)-1;
    end
    variableIdx = unfold.deconv.cols2variableNames(predIDX);
    
    b = beta(:,:,predIDX:predIDX_next);
    
    e = unfold.param(predIDX);
    
    if strcmp(unfold.deconv.variableType{variableIdx},'spline')
        
        splName = cellfun(@(x)x.name,unfold.deconv.splines,'UniformOutput',0);
        splIdx = find(strcmp(splName,unfold.deconv.variableNames{variableIdx}));
        
        spl = unfold.deconv.splines{splIdx};
        % If we found the spline_value given value in the splines, then choose
        % that one
        customSplineValue = strcmp(predNameList,spl.name);
        if any(customSplineValue)
            splValueSelect = predValueSelectList{customSplineValue}{2};
            
        else
            splValueSelect = auto_spacing(cfg,spl.paramValues);
        end
        
        
        % default case, a spline function has been defined
        if isfield(spl,'splinefunction')
            Xspline = spl.splinefunction(splValueSelect,spl.knots);       
            
        elseif spl.knots(1) == spl.knots(2) 
            warning('deprecated spline-function detection,detected default bspline')
            Xspline = default_spline(splValueSelect,spl.knots(3:end-2));
        else
            warning('deprecated spline-function detection,assuming cyclical spline')
            Xspline = cyclical_spline(splValueSelect,spl.knots);
        end
        
        Xspline(:,spl.removedSplineIdx) = [];
        
        
        for c = 1:length(splValueSelect)
            % we have a [channel x time x beta] * [beta x 1] vector product
            % to calculate => loop over channel
            for chan = 1:size(b,1)
                result = permute(b(chan,:,:),[2 3 1])*Xspline(c,:)'; % permute instead of squeeze in case of 1 timepoint 1 channel
                if chan == 1
                    betaNew(chan,:,end+1) =result;
                else
                    betaNew(chan,:,end) =result;
                end
            end
            % the event is already saved in 'e' (unfold.param(predIDX))
            
            eNew = e;
            eNew.value = splValueSelect(c);
            eNew.name = spl.name;
            eNew.type = 'spline_converted';
            
            epochNew(end+1) = eNew;
        end
        
    elseif strcmp(unfold.deconv.variableType{variableIdx},'continuous')
        % we have either a categorical or a continuous predictor here
        customContValue = strcmp(predNameList,unfold.deconv.colnames(predIDX(1)));
        if any(customContValue)
            contValueSelect = predValueSelectList{customContValue}{2};
            
        else
            values = unfold.deconv.X(:,predIDX(1));
            % we do not have access at this point to which 'X' entry belongs to that event.
            % Thus we discard every 0 value and hope a bit, that it does not matter because 0 is the intercept anyway.
            % This is supoptimal and I'm sorry if it creates inconveniences.
            % We would need to introduce a whole new field to carry around
            % to compensate for this.
            warning('auto spacing for continuous variables exlcudes all zeros. Specfiy manually if necessary using ''pred_value''')
            contValueSelect = auto_spacing(cfg,values(values~=0));
        end
        
        for c = contValueSelect
            betaNew(:,:,end+1) = b.*c;
            eNew = e;
            eNew.value = c;
            eNew.type = 'continuous_converted';
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
unfold.param = epochNew;
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
    case 'average'
        contValueSelect = nanmean(predVal);
end
end

