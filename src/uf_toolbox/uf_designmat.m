function [EEG] = uf_designmat(EEG,varargin)
% Generate Designmatrix out of EEG.event structure and a cfg.formula
% Input an EEG event structure and you will get an EEG.unfold.X field with
% the designmatrix.
% If you add multiple eventtypess+formulas as cell-arrays , this function will iteratively
% call itself and combine it to one big designmatrix.
% The designmatrix is not yet ready to do deconvolution, use
% uf_timeexpandDesignmat for this.
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
%                   To define the reference category, have a look at
%                   cfg.categorical {cell} down below. By default matlab
%                   decides on the reference category.
%
%   cfg.eventtypes(cell of strings): cell array of strings, the formula is fit on these
%                  events. make sure that all fields are filled for all events
%                  special-case: Multiple eventtypess (currently in TESTING mode)
%                  You can fit multiple different formulas on different events
%                  concurrently. The specification could be as follows:
%                  {{'A1','A2','A3'},{'B'},{'C'}}.
%                  If more than one formula are specified, we expect you to specify
%                  the eventtypess each formula should be applied to.
%
%   cfg.spline(cell): define a b-spline non-parametric continuous effect. For example:
%                   cfg.spline = {{'speed',15},{'size',10}};
%                   This adds two non-parametric splines, speed with 15
%                   knots and size with 10 knots, thus in total 14+9=23
%                   parameters (we remove one spline due to the Linear
%                   Modeling aspect). Can be specified more conveniently
%                   directly inside the formula
%
%   cfg.categorical(cell-array): default {}, list of which of the EEG.event fields
%                   should be treated as an categorical effect (thus
%                   dummy/effect coded). You can also directly specify what
%                   variables are categorical in the formula.
%                   You can specify the order of the predictors. For
%                   example:
%                   {'predictorA',{'level3','level1','level2'};
%                    'predictorB',{'level2','level1'}}
%                   For predictorA, the level3 is now used as a reference
%                   group. For predictorB the level2 is now used.
%                   The second column of the cell array is optional. E.g.
%                   {'predictorA','predictorB'} will make both predictors
%                   as categorical
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
%     EEG: Returns the EEG structure with the additional fields in EEG.unfold
%
%     * X:          The design matrix
%     * colnames:     For each column of 'X', which predictor it represents
%     * formula:    The original cfg.formula
%     * event:      the cfg.eventtypes
%     * cols2eventtypes:   For each column of 'X' which event it represents
%
%*Example:*
%   A classical 2x2 factorial design with interaction
%|   cfgDesign = [];
%|   cfgDesign.eventtypes = {'fixation'};
%|   cfgDesign.formula = 'y ~ level_predictability*target_fixation + 1';
%|   cfgDesign.categorical = {'level_predictability','target_fixation'};
%|
%|   Second Example
%|   This extends the above example by two cases: A) We add non-parametric
%    splines (n = 10) for the X and Y position of the current fixation. B)
%    We add a second formula for a second event (StimOnset1/2) that only
%    contains a constant (y~1).
%|   cfgDesign.eventtypes = {{'fixation'},{'StimOnset1','StimOnset2'}};
%|   cfgDesign.formula = {'y ~ 1 + level_predictability*target_fixation','y~1'};
%|   cfgDesign.spline = { {{'fixpos_x',10},{'fixpos_y',10}} , {} };
%|   cfgDesign.categorical = {'level_predictability','target_fixation'};
%|
%|   EEG = uf_addDesignmat(EEG,cfgDesign);
%

cfg = finputcheck(varargin,...
    {'categorical',   'cell', [], {};...
    'formula','',[],[];...
    'eventtypes','',[],[];...
    'spline','cell',[],{};...
    'splinespacing','string',{'linear','quantile'},'quantile';
    'codingschema','string',{'effects','reference'},'reference';
    },'mode','ignore');



if(ischar(cfg)); error(cfg);end

if isempty(cfg.eventtypes)
    error('no eventtypes specified even though this is neccesary.')
end

if ~iscell(cfg.eventtypes)
    cfg.eventtypes = {cfg.eventtypes};
end

if iscell(cfg.formula)
    if length(cfg.formula)>1
        display('Multiple events with separate model-formula detected')
        assert(length(cfg.eventtypes) == length(cfg.formula),'You need to specify as many formulas as you have different eventtypes cells')
        % if we are here, eventtypes si something of the likes:
        % {{'A1','A2'},{'B1'},{'C1','C2','C3}} and we need the designmat
        % for each of the cells
        
        
        cfgSingle = cfg;
        for k = 1:length(cfg.eventtypes)
            cfgSingle.eventtypes = cfg.eventtypes{k};
            cfgSingle.formula= cfg.formula{k};
            if ~isempty(cfg.spline)
                cfgSingle.spline = cfg.spline{k};
            end
            
            %do the summersault
            
            EEG2 = uf_designmat(EEG,cfgSingle);
            
            if k == 1
                
                deconvAll = EEG2.unfold;
            else
                %% combine them, rename variables if necessary
                deconvAll = combine_deconvsets(deconvAll,EEG2,k);
            end
            
        end
        EEG.unfold = deconvAll;
        check_rank(EEG.unfold.X)
        %exit the program early
        return
    else
        cfg.formula = cfg.formula{1};
    end
end


% For feedback only
if iscell(cfg.eventtypes) && length(cfg.eventtypes)>1
    eventStr= strjoin(cfg.eventtypes,',');
elseif iscell(cfg.eventtypes)
    eventStr = cfg.eventtypes{1};
else
    eventStr = cfg.eventtypes;
end

fprintf('\nModeling event(s) [%s] using formula: %s \n',eventStr,cfg.formula)
%%
% First of all save the formula
EEG.unfold.formula = {cfg.formula};
%% Regexp the formula
catRegexp = 'cat\((.+?)\)';
splRegexp = 'spl\((.+?)(,.+?)+?\)';

% remove all whitespace
cfg.formula = regexprep(cfg.formula,'[\s]','');

splInteraction= [ regexp(cfg.formula,['2(\*|\:)[\s]*?' splRegexp]) regexp(cfg.formula,[splRegexp '[\s]*?(\*|\:)'])];
if ~isempty(splInteraction)
    error('Spline-Interactions are not supported')
end

cat  = regexp(cfg.formula,catRegexp,'tokens');
spl = regexp(cfg.formula,splRegexp,'tokens');

for s = 1:length(spl)
    if length(spl{s})~=2
        error('error while parsing formula. wrongly defined spline in: %s. Needs to be: spl(your_eventname,10)',cfg.formula)
    end
    spl{s}{2} = str2num(strrep(spl{s}{2},',',''));
    
end%% First check for unique variablenames
cfg.spline = [cfg.spline spl];

% check categorical input and combine
if ~isempty(cfg.categorical)&& size(cfg.categorical,2)>1 && iscell(cfg.categorical{1,2})
    % if one categorical reference exists (a cell array), then all other
    % second entries need to be cell arrays too
    assert(all(cellfun(@(x)iscell(x),cfg.categorical(:,2))),'if you specify one category reference ordering you need to specify all others as well')
    
    categoricalLevels = [cfg.categorical(:,2)];
    
    % the matlab designmatrix function cannot handle {[0],[1]}, we need to
    % transform to {'0' '1'}
    numericix = cellfun(@(x)isnumeric(x{1}),categoricalLevels);
    categoricalLevels(numericix) = cellfun(@(x)cellfun(@(y)num2str(y),x,'UniformOutput',0),categoricalLevels(numericix),'UniformOutput',0);
    
    
    
    categoricalLevelsOrder = [cfg.categorical(:,1)]; %this is needed later to reorder according to the formula
    % we have saved the categoricalLevel information and can now remove it
    % and proceed as planned
    cfg.categorical = [cfg.categorical(:,1)]';
else
    categoricalLevels = []; % matlab should do the ordering
end 


cfg.categorical = unique([cfg.categorical cellfun(@(x)x{1},cat,'UniformOutput',0)]);

% remove the 'cat()' from cat(eventA) so that only 'eventA' remains
f2= regexprep(cfg.formula, catRegexp,'$1');

% remove the spl(splineA) completly
f2= regexprep(f2,['(\+|\*)[\s]*?' splRegexp],'');

% We need to replace the term before the ~ by something that matlab sorts
% after the variables e.g. 'zzz_response'
f2 = regexprep(f2,'.+~','zzz_response~');
% This hack will not be necessary anymore as soon as we have our own parser

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
for k = 1:length(cfg.eventtypes)
    indexmatch = find(strcmpi(cfg.eventtypes(k),{EEG.event.type}));
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
assert(isempty(tf)||all(tf),'Could not find all variables specificed in formula in the data frame generated from the event structure (%s)',sprintf('%s,',F.PredictorNames{~tf}));


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
        % categorical Levels
        %due to a matlab bug, we have to change numeric categorical
        %variables to string categorical variables. Else the reordering of
        %the levels does not work.
        
        if strcmp(t_clean.Properties.VariableNames{p},categoricalLevelsOrder(numericix)) %numericix from line ~210
            is_nan = isnan(t_clean{:,p});
            t_clean.(t_clean.Properties.VariableNames{p}) = arrayfun(@(x)num2str(x),t_clean{:,p},'UniformOutput',0);
            t_clean{is_nan,p} = repmat({''},sum(is_nan),1);
            fprintf('You gave a numeric variable (%s) and specified it as categorical. Due to a matlab-bug we have to convert the variable internally to strings. This should not give further problems \n',t_clean.Properties.VariableNames{p})
        end

        continue
        
    elseif ~iscell(t_clean{:,p}) || ~all(cellfun(@isstr,t_clean{:,p}))
        % We don't have a num (checked before), also not a cell of strings.
        % Error!
        fprintf('non-string events: ')
        fprintf('%i,',find(~cellfun(@isstr,t_clean{:,p})))
        fprintf('\n')
        error('Input Event Values have to be string or numeric. Event:%s was class: %s \n string found in %i out of %i events (sometimes one or a couple of events have NANS instead of strings) ',t.Properties.VariableNames{p},class(t_clean{:,p}),sum(cellfun(@isstr,t_clean{:,p})),size(t_clean,1))
    elseif all(cellfun(@isstr,t_clean{:,p}))
        % If all of them are strings, we need to check that this is a
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
    
    t_clean.zzz_response = zeros(size(t_clean,1),1); % needs to be last to use the designmatrix function
    
    
    is_categorical = ismember(F.VariableNames,cfg.categorical);
    if strcmp(cfg.codingschema,'effects')
        % if we want effects coding, we should remove the mean from the
        % continuous variables as well
        EEG.unfold.effects_mean= nan(1,length(F.VariableNames)-1);% -1 because of the fake-'y'
        for pred = find(~is_categorical(1:end-1))%the last one is the fake-'y'
            predMean = nanmean(t_clean{:,pred});
            EEG.unfold.effects_mean(pred) = predMean;
            t_clean{:,pred} = t_clean{:,pred} - predMean;
        end
        
    end
    
   catlev = cell(size(is_categorical));
    if sum(is_categorical)>0 && ~isempty(categoricalLevels)
        % we need to sort the categoricalLevels to their respective model
        % terms. The user could input the formula in a different order than
        % the cfg..categorical cell-array
        for pIDX = 1:length(categoricalLevelsOrder)
            currCat = categoricalLevelsOrder{pIDX};
            ix = find(strcmp(F.PredictorNames,currCat));
            if isempty(ix)
                error('could not match all cfg.categorical variablenames (%s)to the names in EEG.event',currCat)
            end
            
            %make sure all levels specified occur in the parameter
            dataFrameLevels = unique(t_clean.(currCat))';
            assert(all(cellfun(@(x)any(strcmp(x,dataFrameLevels)),categoricalLevels{pIDX})),'Could not find the levels you specified in cfg.categorical for factor EEG.event.%s ',currCat)
            
            catlev(ix) = categoricalLevels(pIDX);
        end
    end
    categoricalLevels = catlev;

    [X,b,terms,cols2variablenames,colnames] = ...
        classreg.regr.modelutils.designmatrix(t_clean,'Model',F.Terms,...
        'CategoricalVars',is_categorical,'CategoricalLevels',categoricalLevels,...
        'PredictorVars',F.PredictorNames,'ResponseVar','zzz_response',...
        'DummyVarCoding',cfg.codingschema);
    terms = terms';
    %temp = {'continuous','categorical'};
    %predType = {temp{categorical+1}};
end
variablenames = F.VariableNames(1:end-1);
has_intercept = any(strcmp(colnames,'(Intercept)'));

assert(all(cols2variablenames~=0),'error: at least one predictor has only a single value')
if isempty(colnames)
    error('did you specify y~ -1? You need at least a single column in your designmatrix')
end
is_interaction = cellfun(@(x)any(x),strfind(colnames,':'));

if sum(is_interaction)>0
    %check whether main effects were modeled, if not remove them from
    %variablenames
    
    removeList = [];
    
    for int = find(is_interaction)
        predInInteraction = terms(int,:);
        for mainEffectIX = find(predInInteraction)
            % columsn where the maineffect is used
            co = find(terms(:,mainEffectIX));
            if ~any(sum(terms(co,:),2) == 1)
                % we do not have the main effect
                % and remove it
                removeList = [removeList mainEffectIX];
            end
        end
    end
    % Also add the interaction to the variableName List
    variablenames(end+1:(end+sum(is_interaction))) = colnames(is_interaction);
    variablenames(removeList) = [];
    is_categorical(removeList) = [];
end

if has_intercept
    variablenames = ['(Intercept)' variablenames];
end
% cols2variablenames = cols2variablenames - 1; % -1 to remove the intercept which we do not cary explictly in the variablenames (only in the X-colnames)

if any(cols2variablenames)<1
    % I had this problem when adding a categorical predictor that had only
    % one level. I therefore added this check
    error('this error can occur if one of the predictors you specified had only one level/value')
end


%%
% We need to kick out all events that we are not interested in, but keep
% the general event matrix structur the same (e.g. we can still use
% EEG.event(25) and have the entries in EEG.unfold.X(25,:) matching!
% Kicking out ==> putting completly to 0.
X(removeIndex,:) = 0;



EEG.unfold.X = X;



% We want intercept = 0,continuos = 1, categorical = 2, interaction=3 and spline = 4, then add one and
% index the labels

varType = [double(is_categorical(1:end-1)) repmat(2,1,sum(is_interaction))] + 1; %the last categorical one is the fake 'y~',
if has_intercept
    varType = [0 varType];
end

varLabel = {'intercept','continuous','categorical','interaction','spline'};


EEG.unfold.variabletypes = [varLabel(varType+1)];
EEG.unfold.variablenames = variablenames;


EEG.unfold.colnames = colnames;
EEG.unfold.cols2variablenames= cols2variablenames;



EEG.unfold.cols2eventtypes = ones(1,size(X,2)); % This looks odd, but the function is recursively called if multiple events are detected. This number is the fixed in the recursive call to represent the actual cols2eventtypes
EEG.unfold.eventtypes = {cfg.eventtypes};

%% Add the extra defined splines
EEG.unfold.splines = [];
if ~isempty(cfg.spline)
    for s = 1:length(cfg.spline)
        [EEG, ~,nanlist] = uf_designmat_spline(EEG,'name',cfg.spline{s}{1},'nsplines',cfg.spline{s}{2},'paramValues',t{:,cfg.spline{s}{1}},'splinespacing',cfg.splinespacing);
        EEG.unfold.X(nanlist,:) = 0;
    end
end



% Designmat Checks
if any(any(isnan(EEG.unfold.X)))
    warning('NaNs detected in designmat, try to impute them before fitting the model')
    fprintf(['nans found in: ',EEG.unfold.colnames{any(isnan(EEG.unfold.X))}])
    fprintf('\n')
    
end
end

function deconvAll = combine_deconvsets(deconvAll,EEG2,k)
% Function to combine recursive calls
% checks that variables have unique variable names

% Columnnames should be renamend according to variablenames. I.e. the
% following should be:
%   formula: y~condA, y~cat(condA)
%   variablename = 'condA','2_condA'
%   colname = 'condA','2_condA_0','2_condA_1'
% and not
%   colname = 'condA','condA_0','condA_1' 


setA = deconvAll.variablenames;
setA = cellfun(@(x)strsplit(x,':'),setA,'UniformOutput',0);

currnames= EEG2.unfold.variablenames;
colnames = EEG2.unfold.colnames;

for curridx = 1:length(currnames)
    
    currnameext= strsplit(currnames{curridx},':');
    
        % check if any of the extended fieldnames are members of setA
    samenameext = ismember(currnameext,[setA{:}]);
    
    
    
    
    % rename those
    if length(currnameext) == 1 && samenameext
        %maineffect/intercept
       newname = sprintf('%i_%s',k,currnameext{1});
        % which columns belong to this predictor    
        c2var = EEG2.unfold.cols2variablenames;
        
        % find all columns, could be more than variablenames because of
        % dummy coding
        ix = c2var(ismember(c2var,curridx));
        for j = 1:length(ix)
            % find the substring and replace it. Colnames could be longer
            % than variablenames, i.e. 
            % variablename: condA
            % colname:      condA_0, condA_1, condA_2
            colnames{ix(j) + j -1} = strrep(colnames{ix(j) + j - 1},currnameext{1},newname);
        end
        
    elseif any(samenameext)
        % interaction
        newname = [];
        for intidx = 1:length(currnameext)
            if samenameext(intidx)
               newname{intidx} = sprintf('%i_%s',k,currnameext{intidx});
                % colnames
%                 if ismember(j,EEG2.
                  % find all columns, could be more than variablenames because of
                  % dummy coding
                  ix = c2var(ismember(c2var,curridx));
                  for j = ix
                      % find the substring and replace it. Colnames could be longer
                      % than variablenames, i.e.
                      % variablename: condA
                      % colname:      condA_0, condA_1, condA_2
                      colnames{j} = strrep(colnames{j},currnameext{intidx},newname{intidx});
                  end
                  
            else
                newname{intidx} = currnameext{intidx};
            end

        end
       newname = strjoin(newname,':');

    else
        % if we do not need to rename, just take the old name
        newname = currnames{curridx};
    end
    % does not add a ':' for single names

    currnames{curridx} = newname;
    
end
EEG2.unfold.variablenames = currnames;
EEG2.unfold.colnames = colnames;




deconvAll.X(:,(end+1):(end+size(EEG2.unfold.X,2))) = EEG2.unfold.X;
deconvAll.colnames = [deconvAll.colnames,EEG2.unfold.colnames];
deconvAll.formula = [deconvAll.formula EEG2.unfold.formula];
deconvAll.eventtypes = [deconvAll.eventtypes EEG2.unfold.eventtypes];
deconvAll.cols2eventtypes = [deconvAll.cols2eventtypes EEG2.unfold.cols2eventtypes+max(deconvAll.cols2eventtypes)];
deconvAll.variablenames = [deconvAll.variablenames EEG2.unfold.variablenames];
deconvAll.variabletypes = [deconvAll.variabletypes EEG2.unfold.variabletypes];
deconvAll.splines = [deconvAll.splines EEG2.unfold.splines];
%Add the current maximum
EEG2.unfold.cols2variablenames = EEG2.unfold.cols2variablenames+max(deconvAll.cols2variablenames);

% find the ones that where 0 in EEG2.unfold.cols2variablenames
% and reset them to 0 (0 meaning the intercept)
deconvAll.cols2variablenames = [deconvAll.cols2variablenames EEG2.unfold.cols2variablenames];

check_rank(deconvAll.X)
end


function check_rank(X)
    if rank(X)<size(X,2)
        warning('Rank is smaller than matrix size. Do you have two columns that are identical? Other linear dependencies can occur. Check for collinearity')
    end
end