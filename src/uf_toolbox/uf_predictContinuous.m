function output = uf_predictContinuous(ufresult,varargin)
%% Evaluates a continuous/spline parameter at specific values
% This is similar to a predict function, but does not add the marginal of
% the other parameters. For this please make use of uf_addmarginals().
%
% Because model-estimates / parameters are defined for each time-point and
% electrode and can also encompass multiple betas (in the case of spline
% predictors), this becomes non trivial and thus this function.
% Note that this will overwrite the ufresult.beta field
%
%Arguments:
%   cfg.predictAt(cell): One entry per parameter:
%       {{'par1',[10 20 30]},{'par2',[0,1,2]}}.
%       This evaluates parameter 1 at the values 10,20 and 30. Parameter 2
%       at 0, 1 and 2.
%       Default behaviour: evaluates 7 linearly spaced values between the min + max. of the
%       parameterdomain
%   cfg.auto_method(string): 'quantile' (default) or 'linear'.
%       'quantile' - the auto_n values are placed on the quantile of the predictor
%       'linear'   - the auto_n values are placed linearly over the range of the predictor
%       'average'  - only evaluates at the average of the predictor. This
%       is useful if you are interested in the marginal response
%   cfg.auto_n (integer) : default 10; the number of automatically evaluated values
%
%Return:
%   Betas with evaluated betas at specified continuous values.
%
%*Example:*
%    You calculated for a continuous variable "parameterA" a beta of 3.
%    You want to know what the predicted signal of parameterA = [10,20,30] is.
%    You call the function:
%      ufresult = uf_predictContinuous(ufresult,'predictAt',{{'parameterA',[10 20 30]}}
%    The output then would be the respective values 30,60 and 90.


cfg = finputcheck(varargin,...
    {'predictAt','cell',[],{{'',[]}};
    'deconv','integer',[-1,0,1],-1;
    'auto_method','string',{'quantile','linear','average'},'quantile';
    'auto_n','integer',[],10;
    },'mode','error');

if(ischar(cfg)); error(cfg);end

% check whether the user tried to enter EEG.unfold directly into this
% function without running uf_condense first
if ~isfield(ufresult,'param') & isfield(ufresult,'unfold')
    error('\n%s(): You cannot directly enter the unfold output into this function - you have to run uf_condense() first',mfilename)
end

% check if function has been run before
% this will lead to mistakes and errors because some betas will not be the betas we expect anymore - don't allow that!
if any(cellfun(@(x)~isempty(x),strfind({ufresult.param.type},'converted')))
    error('cannot run uf_predictContinuous twice. Run uf_condense again first')
end


beta_dcExists   = isfield(ufresult,'beta')&& isnumeric(ufresult.beta);
beta_nodcExists = isfield(ufresult,'beta_nodc') && isnumeric(ufresult.beta_nodc);

if cfg.deconv == 1
    assert(beta_dcExists,'beta_dc missing or not numeric')
    
    
elseif cfg.deconv== 0
    assert(beta_nodcExists,'beta_nodc missing or not numeric')
    
    
elseif cfg.deconv == -1 % auto detect, recursive call
    
    assert(isfield(ufresult,'beta')|isfield(ufresult,'beta_nodc'),'error: to use autodetect at least the field ufresult.beta  or ufresult.beta_nodc needs to exist')
    fn = fieldnames(ufresult);
    
    if isfield(ufresult,'beta')
        sizeBeta = size(ufresult.beta);
    else
        sizeBeta = size(ufresult.beta_nodc);
    end
    %-------------- Recursive part
    
    for f = fn'
        if strcmp(f,'times')
            continue
        end
        if length(sizeBeta) == length(size(ufresult.(f{1}))) &&  all(sizeBeta == size(ufresult.(f{1})))
            
            % overwrite the "beta" field and use that
            cfg.deconv = 1;
            ufresult_tmp = ufresult;
            ufresult_tmp.beta = ufresult.(f{1});
            output_tmp = uf_predictContinuous(ufresult_tmp,cfg); %------- recursive call
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
predValueSelectList = cfg.predictAt;
predNameList = cellfun(@(x)x{1},predValueSelectList,'UniformOutput',0);
[~,paramList] = uf_getSplineidx(ufresult);
if cfg.deconv == 1
    beta = ufresult.beta;
elseif cfg.deconv== 0
    
    beta = ufresult.beta_nodc;
end

betaNew = [];
epochNew = ufresult.param(1); %needs to be removed


for currPred = 1:length(paramList)
    predIDX = paramList(currPred);
    
    % if we are at the last predictor, this one goes to the end of the
    % designmatrix
    if currPred == length(paramList)
        predIDX_next = size(ufresult.unfold.X,2);
    else
        %else it goes to the next predictor
        predIDX_next = paramList(currPred+1)-1;
    end
    variableIdx = ufresult.unfold.cols2variablenames(predIDX);
    
    b = beta(:,:,predIDX:predIDX_next);
    
    e = ufresult.param(predIDX);
    
    if strcmp(ufresult.unfold.variabletypes{variableIdx},'spline')
        
        splName = cellfun(@(x)x.name,ufresult.unfold.splines,'UniformOutput',0);
        splIdx = find(strcmp(splName,ufresult.unfold.variablenames{variableIdx}));
        
        spl = ufresult.unfold.splines{splIdx};
        % If we found the spline_value given value in the splines, then choose
        % that one
        customSplineValue = strcmp(predNameList,spl.name);
        if any(customSplineValue)
            if size(spl.paramValues,1) == 2
                splValueSelect = predValueSelectList{customSplineValue}(2:3);
            else
                splValueSelect = predValueSelectList{customSplineValue}{2};
            end
        else
            splValueSelect = auto_spacing(cfg,spl.paramValues);
        end
        
        
        % default case, a spline function has been defined
        if isfield(spl,'splinefunction')
            if size(spl.paramValues,1) == 2
                auxsiz = size(splValueSelect{2},2);
                splValueSelect{2}  = reshape(repmat(splValueSelect{2},length(splValueSelect{1}),1),1,length(splValueSelect{1})* size(splValueSelect{2},2));
                splValueSelect{1}  = repmat(splValueSelect{1},1,auxsiz);
                Xspline1 =  spl.splinefunction(splValueSelect{1},spl.knots(1,:));
                Xspline2 =  spl.splinefunction(splValueSelect{2},spl.knots(2,:));
                Xspline  = [];
                for iD = 1:size(Xspline1,1)
                    aux = Xspline1(iD,:)'*Xspline2(iD,:);
                    Xspline(iD,:) = aux(:)';
                end
            else
                Xspline = spl.splinefunction(splValueSelect,spl.knots);
            end
        elseif spl.knots(1) == spl.knots(2)
            warning('deprecated spline-function detection,detected default bspline')
            Xspline = default_spline(splValueSelect,spl.knots(3:end-2));
        else
            warning('deprecated spline-function detection,assuming cyclical spline')
            Xspline = cyclical_spline(splValueSelect,spl.knots);
        end
        
        Xspline(:,spl.removedSplineIdx) = [];
        
        for c = 1:size(Xspline,1)
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
            % the event is already saved in 'e' (ufresult.param(predIDX))
            
            eNew = e;
            if iscell(splValueSelect)
                eNew.value = [splValueSelect{1}(c),splValueSelect{2}(c)];
            else
                eNew.value = splValueSelect(c);
            end
            eNew.name = spl.name;
            eNew.type = 'spline_converted';
            
            epochNew(end+1) = eNew;
        end
        
    elseif strcmp(ufresult.unfold.variabletypes{variableIdx},'continuous')
        
        % we have either a categorical or a continuous predictor here
        customContValue = strcmp(predNameList,ufresult.unfold.colnames(predIDX(1)));
        if any(customContValue)
            contValueSelect = predValueSelectList{customContValue}{2};
        else
            values = ufresult.unfold.X(:,predIDX(1));
            % we do not have access at this point to which 'X' entry belongs to that event.
            % Thus we discard every 0 value and hope a bit, that it does not matter because 0 is the intercept anyway.
            % This is supoptimal and I'm sorry if it creates inconveniences.
            % We would need to introduce a whole new field to carry around
            % to compensate for this.
            warning(sprintf('The predictor ''%s'' is evaluated at its mean or at (auto-spaced) percentiles. Note that the computation of the mean/percentiles for continuous variables excludes all cases where the predictor value is exactly zero. If necessary, specify the values at which to evaluate the predictor manually using ''predictAt''',ufresult.unfold.variablenames{variableIdx}))
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
    ufresult.beta = betaNew;
else
    ufresult.beta_nodc = betaNew;
end
ufresult.param = epochNew;
output = ufresult;

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
        warning('moving min/max inside by 0.05 of total range in order to reduce bad extreme-estimates');
        contValueSelect = linspace(contmin+0.05*ran,contmax-0.05*ran,cfg.auto_n);
    case 'quantile'
        %contValueSelect = quantile(predVal,linspace(0,1,cfg.auto_n));
        %contValueSelect = quantile(predVal,linspace(1/cfg.auto_n,1-1/cfg.auto_n,cfg.auto_n-1)); % VERSION OLAF
        contValueSelect = quantile(predVal,linspace( 1/(cfg.auto_n+1), 1-1/(cfg.auto_n+1), cfg.auto_n)) % update Olaf
    case 'average'
        contValueSelect = nanmean(predVal);
end
end

