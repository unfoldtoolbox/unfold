function [varargout] = dc_plotParam(unfold,varargin)
% Plots time vs. Voltage in separate plots for each predictor, where there
% are multiple lines for each predictor
%
% 'unfold' needs to have the 'unfold' structure, the output from
% "dc_beta2unfold"
%
% Uses the 'gramm'-toolbox for plotting
%
%Arguments:
%    'channel' (integer): Which channel to plot
%
%    'pred_value' (cell): a cell of cell arrays, e.g. {{'parName',linspace(0,10,5)},{'parname2',1:5}}
%       This splits up the parName-predictor into 5 bins from
%       0 to 10, so 5 lines would be plotted. Default are 7 lines
%       from min to max
%
%    'deconv' ([-1 0 1]):default: -1; whether to plot unfold.beta (1) or
%       unfold.beta_nodc(0) or everything/autodetect (-1). Autodetect would
%       also detect same-shaped other predictors. If e.g. you want to compare
%       multiple runs from different algorithms or similar
%
%    'add_intercept' (boolean): Add the intercept/constant to each subplot.
%      This will give ERP-plots that are commonly used. Without add_intercepts the factors (if they are categorical) could be interpretet as difference or sometimes main effect plots (if effects-coding is used)
%
%    'baseline' (2 integers): default none; Performs a baseline corrections on the interval (in seconds = unfold.times units) given.
%
%    'include_intercept' (boolean) : default 0; useful with "add_intercept", will add the constant/intercept to each subplot
%
%    'plotSeparate' ('all','event','none'): Each predictor will be
%      plotted in a separate figure ('all'), plotted in an event-specific
%      figure ('event') or all subplots are in the same figure ('none', default)
%
%    'plotParam' (string/cell of strings): Defines which parameters are to be plotted
%
%    'sameyaxis' ('all','row','independent'): Force the same y-axis
%      (default 'all')
%
%    'gramm': (gramm-object) plots the current data ontop of the last gramm-object. This is
%      useful to plot multiple subjects in a single figure.
%
%    'figure' (boolean): Generate a new figure? (default 1)
%
%Returns:
%   allAxesInFigure: All 'subplot' axes that were generated
%
%Example
% dc_plotParam(unfold,'channel',1)

% parse inputs
cfg = finputcheck(varargin,...
    {'pred_value','cell',[],{{'',[]}};
    'deconv','integer',[-1 0 1],-1;
    'channel','',[],[];
    'add_intercept','boolean',[],0;
    'add_marginal', 'boolean', [],0;
    'include_intercept','boolean',[],0;
    'plotParam','',[],'';
    'plotSeparate','string',{'all','event','none'},'none';
    'sameyaxis','string',{'all','row','independent'},'all';
    'baseline','real',[min(unfold.times) max(unfold.times)],[];...
    'gramm','','',[];
    'figure','boolean',[0 1],1;
    },'mode','ignore');
if(ischar(cfg)); error(cfg);end

if ischar(cfg.channel)
   assert(~isempty(unfold.chanlocs),'unfold.chanlocs is empty, it is necessary to be non-empty if you want to specify a channel by string, use numbers instead or populate unfold.chanlocs') 
   cfg.channel = find(strcmp({unfold.chanlocs.labels},cfg.channel));
end
assert(~isempty(cfg.channel)&& cfg.channel>0 &&length(cfg.channel) == 1,'error please select a single channel to plot')

assert(~(cfg.add_marginal&&cfg.add_intercept),'cannot add average AND intercept (the former contains the latter')

% Find out whether we want beta_dc, beta_nodc and if there are other fields
% that have the same size that we should plot as columns.
betaSetName = dc_unfoldbetaSetname(unfold,varargin{:});

if ~isempty(cfg.pred_value{1}{1})
    fprintf('Evaluating parameters at auto or specified values');
    unfold = dc_getParam(unfold,cfg);
end
%% Prepare data
% select parameters to plot, or else plot all available

if isempty(cfg.plotParam)
    display('plotting all parameters')
    paramList = {unfold.param.name};
    paramIdx = 1:length(unfold.param);
else
    display('plotting selected parameters')
    paramIdx = [];
    if isstr(cfg.plotParam)
        % if only a single parameter is requested
        cfg.plotParam = {cfg.plotParam};
    end
    paramList = {};
    for c = cfg.plotParam
        hits = find(strcmp(c, {unfold.param.name}));
        paramIdx = [paramIdx hits];
        paramList = [paramList {unfold.param(hits).name}];
    end
    
end

% add the variable names for the plot
%name + event are later used to split up columns

name = paramList;
event = cellfun(@(y)strjoin_custom(y,'+'),{unfold.param(paramIdx).event},'UniformOutput',0);

% value, this is necessary for continuous + spline
value = [unfold.param(paramIdx).value];
value(isnan(value)) = 0; %categorical variables are nan, we need to convert


if cfg.add_marginal
    % necessary to add the average of the splines, i.e. simulate a marginal plot
    unfold = dc_addmarginal(unfold,'channel',cfg.channel);
end

%the linestyle, it is used when the intercept is added to differentiate
%between intercept and actual factor
plotData = []; plotLinestyle = {}; plotColLabel = [];plotName = [];plotEvent=[];plotValue = [];
for bName =betaSetName
    data = permute(unfold.(bName{1})(cfg.channel,:,paramIdx),[2 3 1]); % squeeze transposes, if paramIDX and channel is 1
    
    if ~isempty(cfg.baseline)
        fprintf('performing baseline correction \n')
        data = bsxfun(@minus,data,mean(data((unfold.times>=cfg.baseline(1))& (unfold.times<cfg.baseline(2)),:),1));
    end
    if length(size(unfold.(bName{1}))) == 2
        data = data';
    end
    
    
    
    plotLinestyle = [plotLinestyle repmat({'Predictor'},1,size(data,2))];
    plotColLabel = [plotColLabel repmat(bName,1,size(data,2))];
    plotName = [plotName paramList];
    plotEvent = [plotEvent event];
    plotValue = [plotValue value];
    
    if cfg.add_intercept || cfg.include_intercept
        
        for e  = unique(event)
            %is there an intercept?
            eventIdx =  strcmp(event,e);
            
            interceptIdx = cellfun(@(x)~isempty(x),strfind({unfold.param(paramIdx).name},'(Intercept)'));
            if sum(interceptIdx) == 0
                warning('no intercept found, did you select a parameter but not the intercepts?')
                continue
            end
            if sum(eventIdx&interceptIdx)>1
                error('multiple intercepts per event are impossible')
            end
            %extract it
            betaIntercept = data(:,eventIdx&interceptIdx);
            
            
            addToIdx = find(eventIdx&~interceptIdx);
            
            if cfg.add_intercept
                %add it
                data(:,addToIdx) = bsxfun(@plus,data(:,addToIdx),betaIntercept);
            end
            if cfg.include_intercept
                % This is mostly useful for categorical factors
                % we need another entry for each subplot, thus we copy
                % everything over
                plotLinestyle = [plotLinestyle repmat({'intercept'},1,length(addToIdx))];
                plotName = [plotName name(addToIdx)];
                plotEvent = [plotEvent event(addToIdx)];
                plotValue = [plotValue zeros(1,length(addToIdx))];
                plotColLabel= [plotColLabel repmat(bName,1,length(addToIdx))];
                
                % we add the data in the end
                data(:,(end+1):(end+length(addToIdx))) =  repmat(betaIntercept,1,length(addToIdx));
                
                
            end
        end
    end
    
    plotData = cat(2,plotData,data);
    
end

%% Make figure


if cfg.figure && isempty(cfg.gramm)
    
    if isfield(unfold,'chanlocs') && ~isempty(unfold.chanlocs)
        channame = unfold.chanlocs(cfg.channel).labels;
    else
        channame = num2str(cfg.channel);
    end
    figure('name',sprintf('unfold-toolbox channel %s',channame))
end
clear g



% Combine name + event incase we plot everything in one plot
plotLabel = cellfun(@(x,y)[x '@' y],plotName,plotEvent,'UniformOutput',0);
rowName = plotLabel;

% Should all axis be the same?
switch cfg.sameyaxis
    case 'all'
        facet_scale = 'fixed';
    case 'row'
        facet_scale = 'free_y';
    case 'independent'
        facet_scale = 'independent';
end

% Generate gramm, update old one if it exists
if isempty(cfg.gramm)
    g = gramm('x',unfold.times,'y',plotData','group',plotLabel,'color',plotValue,'linestyle',plotLinestyle);
else
    g = cfg.gramm;
    update(g,'x',unfold.times,'y',plotData','group',plotLabel,'color',plotValue,'linestyle',plotLinestyle);
end

% which variables should be plotted in new figures?
switch cfg.plotSeparate
    case 'all'
        g.fig(rowName);
        rowName = plotName;
    case 'event'
        g.fig(plotEvent);
        rowName = plotName;
    case 'none'
        rowName = rowName;
        
end

% the actuall plotting things
facet_grid(g,plotColLabel,rowName,'scale',facet_scale);
geom_line(g);

%give us continuous colors!
if length(unique(plotValue))>3
    set_continuous_color(g,'colormap','viridis');
end
if isempty(cfg.gramm)
    geom_vline(g,'xintercept',0);
    geom_hline(g,'yintercep',0);
end
set_names(g,'color','','x','Time [s]','y','ERP [yV]','linestyle','','row','','column','');

draw(g);

if nargout == 1
    varargout{1} = g;
end
