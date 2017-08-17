function [EEG] = dc_designmat(EEG,varargin)
% Generate Designmatrix out of EEG.event structure and a cfg.formula
% Input an EEG event structure and you will get an EEG.deconv.X field with
% the designmatrix.
% If you add multiple eventtypes+formulas as cell-arrays , this function will iteratively
% call itself and combine it to one big designmatrix.
% The designmatrix is not yet ready to do deconvolution, use
% dc_timeexpandDesignmat for this.
%
%Arguments:
%
%   cfg.formula(string):     Formula in the wilkinson format. In addition
%                   one can specify 'cat(X)' so that X is interpreted as a
%                   categorical variable and therefore dummy/effect-coded.
%                   Also using spl(Y,5) defines a "non-linear" predictor
%                   using 5 b-cubic spline basis functions. The more
%                   splines one uses, the higher the risk of overfitting
%                   but also of course more flexible relations can be
%                   fitted.
%
%                   Example with multiple formulas: {'y~A+spl(B,5)', 'y~x+cat(y)','y~1'}
%
%                   Example with more complex formula:
%                   {'y ~ stimulus_type + color * size + stimulus_type:color}'
%                   This formula would add the following main effects:
%                   "Stimulus Type, color, size"
%                   and the following interactions:
%                   "stimType:color, color:size"
%
%   cfg.eventtype(cell of strings): cell array of strings, the formula is fit on these
%                  events. make sure that all fields are filled for all events
%                  special-case: Multiple eventtypes (currently in TESTING mode)
%                  You can fit multiple different formulas on different events
%                  concurrently. The specification could be as follows:
%                  {{'A1','A2','A3'},{'B'},{'C'}}.
%                  If more than one formula are specified, we expect you to specify
%                  the eventtypes each formula should be applied to.
%
%   cfg.spline(cell): define a b-spline non-parametric continuous effect. For example:
%                   cfg.spline = {{'speed',15},{'size',10}};
%                   This adds two non-parametric splines, speed with 15
%                   knots and size with 10 knots, thus in total 14+9=23
%                   parameters (we remove one spline due to the Linear
%                   Modeling aspect). Can be specified more conveniently
%                   directly inside the formula
%
%   cfg.categorical(array): default [], which of the EEG.event fields
%                   should be treated as an categorical effect (thus
%                   dummy/effect coded). You can also directly specify what
%                   variables are categorical in the formula.
%
%   cfg.splinespacing (string): defines how the knots of the splines should be
%                  placed. Possible values:
%                  'linear' : linear spacing with boundary splines at the respective min/max
%                  'log': logarithmic increasing spacing
%                  'logreverse': log decreasing spacing
%                  'quantiles' (default): heuristic spacing at the quantiles
%
%   cfg.codingschema(string): default: 'references', could be 'effects', this is
%                   relevant if you define categorical input variables.
%                   Reference coding is also known as treatment coding
%
%Returns:
%     EEG: Returns the EEG structure with the additional fields in EEG.deconv
%
%     * X:          The design matrix
%     * colnames:     For each column of 'X', which predictor it represents
%     * formula:    The original cfg.formula
%     * event:      the cfg.eventtype
%     * col2eventtype:   For each column of 'X' which event it represents
%
%*Example:*
%   A classical 2x2 factorial design with interaction
%|   cfgDesign = [];
%|   cfgDesign.eventtype = {'fixation'};
%|   cfgDesign.formula = 'y ~ level_predictability*target_fixation + 1';
%|   cfgDesign.categorical = {'level_predictability','target_fixation'};
%|
%|   Second Example
%|   This extends the above example by two cases: A) We add non-parametric
%    splines (n = 10) for the X and Y position of the current fixation. B)
%    We add a second formula for a second event (StimOnset1/2) that only
%    contains a constant (y~1).
%|   cfgDesign.eventtype = {{'fixation'},{'StimOnset1','StimOnset2'}};
%|   cfgDesign.formula = {'y ~ 1 + level_predictability*target_fixation','y~1'};
%|   cfgDesign.spline = { {{'fixpos_x',10},{'fixpos_y',10}} , {} };
%|   cfgDesign.categorical = {'level_predictability','target_fixation'};
%|
%|   EEG = dc_addDesignmat(EEG,cfgDesign);
%

cfg = finputcheck(varargin,...
    {'categorical',   'cell', [], {};...
    'formula','',[],[];...
    'eventtype','',[],[];...
    'spline','cell',[],{};...
    'splinespacing','string',{'linear','log','quantile','logreverse'},'quantile';
    'codingschema','string',{'effects','reference'},'reference';
    },'mode','ignore');



if(ischar(cfg)); error(cfg);end

if isempty(cfg.eventtype)
    error('no eventtype specified even though this is neccesary.')
end

if ~iscell(cfg.eventtype)
    cfg.eventtype = {cfg.eventtype};
end

if iscell(cfg.formula)
    if length(cfg.formula)>1
        display('Multiple events with separate model-formula detected')
        assert(length(cfg.eventtype) == length(cfg.formula),'You need to specify as many formulas as you have different eventtype cells')
        % if we are here, eventtype si something of the likes:
        % {{'A1','A2'},{'B1'},{'C1','C2','C3}} and we need the designmat
        % for each of the cells
        
        
        cfgSingle = cfg;
        for k = 1:length(cfg.eventtype)
            cfgSingle.eventtype = cfg.eventtype{k};
            cfgSingle.formula= cfg.formula{k};
            if ~isempty(cfg.spline)
                cfgSingle.spline = cfg.spline{k};
            end
            
            %do the summersault
            
            EEG2 = dc_designmat(EEG,cfgSingle);
            if k == 1
                
                deconvAll = EEG2.deconv;
            else
                deconvAll.X(:,(end+1):(end+size(EEG2.deconv.X,2))) = EEG2.deconv.X;
                deconvAll.colnames = [deconvAll.colnames,EEG2.deconv.colnames];
                deconvAll.formula = [deconvAll.formula EEG2.deconv.formula];
                deconvAll.eventtype = [deconvAll.eventtype EEG2.deconv.eventtype];
                deconvAll.col2eventtype = [deconvAll.col2eventtype EEG2.deconv.col2eventtype+max(deconvAll.col2eventtype)];
                deconvAll.variableNames = [deconvAll.variableNames EEG2.deconv.variableNames];
                deconvAll.variableType = [deconvAll.variableType EEG2.deconv.variableType];
                deconvAll.predictorSplines = [deconvAll.predictorSplines EEG2.deconv.predictorSplines];
                %Add the current maximum
                EEG2.deconv.cols2variableNames = EEG2.deconv.cols2variableNames+max(deconvAll.cols2variableNames);
                
                % find the ones that where 0 in EEG2.deconv.cols2variableNames
                % and reset them to 0 (0 meaning the intercept)
                ix = EEG2.deconv.cols2variableNames==max(deconvAll.cols2variableNames);
                EEG2.deconv.cols2variableNames(ix) = 0;
                deconvAll.cols2variableNames = [deconvAll.cols2variableNames EEG2.deconv.cols2variableNames];
            end
            
        end
        EEG.deconv = deconvAll;
        %exit the program early
        return
    else
        cfg.formula = cfg.formula{1};
    end
end


% For feedback only
if iscell(cfg.eventtype) && length(cfg.eventtype)>1
    eventStr= strjoin(cfg.eventtype,',');
elseif iscell(cfg.eventtype)
    eventStr = cfg.eventtype{1};
else
    eventStr = cfg.eventtype;
end

fprintf('Modeling event(s) [%s] using formula: %s \n',eventStr,cfg.formula)

%% Regexp the formula
catRegexp = 'cat\((.+?)\)';
splRegexp = 'spl\((.+?)(,.+?)+?\)';


splInteraction= [ regexp(cfg.formula,['2(\*|\:)[\s]*?' splRegexp]) regexp(cfg.formula,[splRegexp '[\s]*?(\*|\:)'])];
if ~isempty(splInteraction)
    error('Spline-Interactions are not supported')
end




splRegexp = 'spl\((.+?)(,.+?)+?\)';
cat  = regexp(cfg.formula,catRegexp,'tokens');
spl = regexp(cfg.formula,splRegexp,'tokens');

for s = 1:length(spl)
    spl{s}{2} = str2num(strrep(spl{s}{2},',',''));
end
cfg.spline = [cfg.spline spl];

cfg.categorical = [cfg.categorical cellfun(@(x)x{1},cat,'UniformOutput',0)];
f2= regexprep(cfg.formula, catRegexp,'$1');
f2= regexprep(f2,['(\+|\*)[\s]*?' splRegexp],'');
cfg.formula = f2;
%%

% Parse the formula
F = classreg.regr.LinearFormula(cfg.formula);

% init the table for the desigmnat function

event = EEG.event;

% go through all events and replace [] with nan
% else the struct2table will make cells out of numerics
fn = fieldnames(event);
for e= 1:length(event)
    for f = fn'
        if isempty(event(e).(f{1}))
            event(e).(f{1}) = nan;
        end
    end
end
t = struct2table(event);




% find and remove all event types that are not needed
% deletes **** rows ****
indexList = [];
for k = 1:length(cfg.eventtype)
    indexmatch = find(strcmpi(cfg.eventtype(k),{EEG.event.type}));
    if isempty(indexList)
        indexList = indexmatch;
    else
        indexList = union(indexList,indexmatch);
    end
    %
end
removeIndex = setdiff(1:size(t,1),indexList);

% we want only the designmatrix of triggers we want to model
for p = 1:size(t,2)
    %this awkward construction is entirely matlabs fault ;)
    try
        t{removeIndex,p} =nan;
        t(:,p) = standardizeMissing(t(:,p),nan);
        
    catch
        t{removeIndex,p} ={''};
    end
end

%

% Check everything necessary is there
tf = ismember(F.PredictorNames,t.Properties.VariableNames);
assert(all(tf),'Could not find all variables specificed in formula in the data frame generated from the event structure (%s)',sprintf('%s,',F.PredictorNames{~tf}));


% Remove every variable we dont need
% deletes **** columns****
[tf] = find(ismember(t.Properties.VariableNames,F.PredictorNames));

% bring them in the right order
reorder = [];
for k = 1:length(F.PredictorNames)
    reorder = [reorder find(strcmp(t.Properties.VariableNames(tf),F.PredictorNames(k)))];
end

tf = tf(reorder);

t_clean = t(:,tf);

%Check if all input is strings or numeric
% numeric input is easy

t_isnum = varfun(@isnumeric,t_clean,'output','uniform');
for p = 1:size(t_clean,2)
    if t_isnum(p)
        
        continue
        
    elseif ~iscell(t_clean{:,p}) || ~all(cellfun(@isstr,t_clean{:,p}))
        % We don't have a num (checked before), also not a cell of strings.
        % Error!
        error('Input Event Values have to be string or numeric. Event:%s was class: %s',t.Properties.VariableNames{p},class(t_clean{:,p}))
    elseif all(cellfun(@isstr,t_clean{:,p}))
        % If all of them arr strings, we need to check that this is a
        % categorical variable
        currVar = t_clean.Properties.VariableNames{p};
        if ~ismember(currVar,cfg.categorical)
            warning('The specified variable: "%s", was detected as a string but not marked as a categorical variable. We model it as categorical, be sure this is what you want.',currVar)
            cfg.categorical{end+1} = currVar;
        end
    end
end


% add a dummy response variable
if isempty(F.Terms)% no variables selected, maybe only splines?
    fprintf('No parametric effects selected, Only splines are generated\n')
    X = [];
    colnames = [];
    %predType = [];
else
    
    t_clean.y = zeros(size(t_clean,1),1); % needs to be last to use the designmatrix function
    
    
    categorical = ismember(F.VariableNames,cfg.categorical);
    if strcmp(cfg.codingschema,'effects')
        % if we want effects coding, we should remove the mean from the
        % continuous variables as well
        EEG.deconv.effects_mean= nan(1,length(F.VariableNames)-1);% -1 because of the fake-'y'
        for pred = find(~categorical(1:end-1))%the last one is the fake-'y'
            predMean = nanmean(t_clean{:,pred});
            EEG.deconv.effects_mean(pred) = predMean;
            t_clean{:,pred} = t_clean{:,pred} - predMean;
        end
        
    end
    
    [X,terms,a,cols2variableNames,colnames] = ...
        classreg.regr.modelutils.designmatrix(t_clean,'Model',F.Terms,...
        'CategoricalVars',categorical,...
        'PredictorVars',F.PredictorNames,'ResponseVar','y',...
        'DummyVarCoding',cfg.codingschema);
    %temp = {'continuous','categorical'};
    %predType = {temp{categorical+1}};
end
variableNames = F.VariableNames(1:end-1);
has_intercept = any(strcmp(colnames,'(Intercept)'));
if has_intercept
    variableNames = ['(Intercept)' variableNames];
end
if isempty(colnames)
    error('did you specify y~ -1? You need at least a single column in your designmatrix')
end
is_interaction = cellfun(@(x)any(x),strfind(colnames,':'));

if sum(is_interaction)>0
    variableNames(end+1:(end+sum(is_interaction))) = colnames(is_interaction);
end
% cols2variableNames = cols2variableNames - 1; % -1 to remove the intercept which we do not cary explictly in the variableNames (only in the X-colnames)

if any(cols2variableNames)<1
    % I had this problem when adding a categorical predictor that had only
    % one level. I therefore added this check
    error('this error can occur if one of the predictors you specified had only one level/value')
end
%% Add the extra defined splines. In the end I want something like y ~ x1 + x2*x3 + s(x4) to be parsed correctly... but this is more difficult to implement
EEG.deconv.predictorSplines = [];
if ~isempty(cfg.spline)
    
    for s = 1:length(cfg.spline)
        spl = [];
        spl.paramValues = t{:,cfg.spline{s}{1}};
        spl.nSplines = cfg.spline{s}{2};
        %spl.range = range(spl.paramValues);
        spl.min = min(spl.paramValues);
        spl.max= max(spl.paramValues);
        spl.name = [cfg.spline{s}{1}];
        % We find the smallest increase to get an idea of the resolution we
        % need to define the splines
        
        nanlist = [];
        
        spl.spline2val = linspace(spl.min,spl.max,100*spl.nSplines);
        %spl.spline2val = spl.min:min(b):spl.max;
        for k = 1:length(spl.paramValues)
            if isnan(spl.paramValues(k)) && ~isempty(t.type{k})
                
                warning('nan found @ %s, event:%i \n currently no support for nan-value in splines. setting it to 0.',cfg.spline{s}{1},k)
                
                nanlist = [nanlist k];
                spl.spline2val_idx(k)=1; % this is only temporary, we will remove it after indexing again
                continue
            end
            [~,spl.spline2val_idx(k)] = (min(abs(spl.spline2val- spl.paramValues(k))));
        end
        knots = [];
        
        
        switch cfg.splinespacing
            case 'linear'
                knots = linspace(1,length(spl.spline2val),spl.nSplines-2);
            case 'log'
                knots = logspace(log10(1),log10(length(spl.spline2val)),spl.nSplines-2);
            case 'logreverse'
                knots = logspace(log10(1),log10(length(spl.spline2val)),spl.nSplines-2);
                knots = sort(length(spl.spline2val) - (knots-1));
            case 'quantile'
                % here we are not working on 0 to X (e.g. X could be 1500 = nSplines*100) but we calculate the
                % quantiles on the predictor data
                knots = quantile(spl.paramValues,linspace(0,1,spl.nSplines-2));
                % now we have to find where the quantiles are in the 0:X
                % space.
                for k = 1:length(knots)
                    [~,knots(k)] = (min(abs(spl.spline2val- knots(k))));
                end
        end
        % we add 3 knots (because cubic, 4th order splines) in the
        % beginning and the end.
        knots = [repmat(knots(1),1,3) knots repmat(knots(end),1,3)];
        spl.basis = Bernstein(1:length(spl.spline2val),knots,[],4); % 4 is the order
        spl.basis(end,end) = 1; % there seems to be a bug in the above function. The last entry should be 1 and not 0
        
        %%
        if ~strcmp(cfg.codingschema,'effects')
            warning('For splines only effects coding is implemented. The intercept will represent the value at the spline most closely to the mean')
        end
        [~,killThisSpline] = max(spl.basis(get_min(nanmean(spl.paramValues),spl.spline2val),:)); % we remove the spline that is closed to the mean value, thus it will kind of be like an effects coding?
        fprintf('The intercept for the effect %s has its peak at %f\n',spl.name,nanmean(spl.paramValues))
        %         else
        %             killThisSpline = 1;
        %         end
        spl.basis(:,killThisSpline) = [];
        
        spl.X= spl.basis(spl.spline2val_idx,:);
        spl.X(nanlist,:) = 0; % remove the nans
        X(nanlist,:) = 0; % and remove them from the designmatrix (see warning above)
        
        spl.nSplines = size(spl.X,2); % knots might change due to the create_splines_linspace internal function
        X = [X spl.X];
        
        
        
        rawColnames = repmat({spl.name},1,spl.nSplines);
        % find the max of each spline (the value that best represents this
        % predictor-term
        [~,I] = max(spl.basis,[],1);
        maxSplVal = spl.spline2val(I);
        
        % round to two significant digits:
        tmpSplVal = 2-floor(log10(abs(maxSplVal)));
        tmpSplVal(tmpSplVal<0) = 0;
        spl.colnames = cellfun(@(x,signPoint,y)sprintf('%s_%.*f',x,signPoint,y),rawColnames,num2cell(tmpSplVal),num2cell(maxSplVal),'UniformOutput', 0);
        colnames = [colnames  spl.colnames]; % TODO change the name so that each column depicts the mean value of the spline in addition
        variableNames = [variableNames {spl.name}];
        
        
        cols2variableNames = [cols2variableNames repmat(length(variableNames),1, spl.nSplines)];
        %predType = [predType repmat({'spline'},1,spl.nSplines)];
        EEG.deconv.predictorSplines{end+1} = spl;
    end
end

% We need to kick out all events that we are not interested in, but keep
% the general event matrix structur the same (e.g. we can still use
% EEG.event(25) and have the entries in EEG.deconv.X(25,:) matching!
% Kicking out ==> putting completly to 0.
X(removeIndex,:) = 0;

EEG.deconv.formula = {F};

EEG.deconv.X = X;



% We want intercept = 0,continuos = 1, categorical = 2, interaction=3 and spline = 4, then add one and
% index the labels

varType = [double(categorical(1:end-1)) repmat(2,1,sum(is_interaction)) repmat(3,1,length(cfg.spline)) ] + 1; %the last categorical one is the fake 'y~', 
if has_intercept
    varType = [0 varType];
end

varLabel = {'intercept','continuous','categorical','interaction','spline'};


EEG.deconv.variableType = [varLabel(varType+1)];
EEG.deconv.variableNames = variableNames;


EEG.deconv.colnames = colnames;
EEG.deconv.cols2variableNames= cols2variableNames;



EEG.deconv.col2eventtype = ones(1,size(X,2)); % This looks odd, but the function is recursively called if multiple events are detected. This number is the fixed in the recursive call to represent the actual col2eventtype
EEG.deconv.eventtype = {cfg.eventtype};


%EEG.deconv.predType = predType;

% Designmat Checks
if any(any(isnan(EEG.deconv.X)))
    warning('NaNs detected in designmat, be sure to impute them before fitting the model')
    warning(['found in: ',EEG.deconv.colnames{any(isnan(EEG.deconv.X))}])
    
end
