function [varargout] = uf_plotParamTopo(ufresult,varargin)
%Generates topoplots over time. each row is a predictor or predictor-value-combination
%
%Arguments:
% 'plotParam' : cell array of parameters to be plotted, if empty plots all
% 'predictAt': only used for continuous / spline predictors. 
%               a cell of cell arrays, e.g. {{'parName',linspace(0,10,5)},{'parname2',1:5}}
%               This would split up the parName-predictor into 5 bins from
%               0 to 10, so 5 rows of topoplots would be plotted. Default are 7 lines
%               from min to max (see code for more info)
%
% 'add_intercept' : (not yet supported) whether to add the intercept to each curve (thus not
%               the 'pure' effect of a independent variable is plotted but
%               more of an ERP-like plot is generated
%
% 'n_topos' :(15) number of topographies to plot
%
% 'channel' : plot only a subset of channels
%
% 'baseline' (2 integers): default none; Performs a baseline corrections on the interval (in seconds = ufresult.times units) given.
%
% 'caxis' ('same',default:[]) if 'same', generates the same coloraxis based
%       on the 95% percentile of the selected beta-values. can be
%       customized to whichever caxis e.g. [-3 5]
% 'betaSetName'   : Default 'beta'. Can be any field of the ufresult-struct
%
% 'figure'    : plot in new figure (1) or old (0), default: (1)
%
%Returns:
%   structure of all plotting axes.
%
%*Examples:*
%  uf_plotParamTopo(EEG,'plotParam',{'FactorX','FactorC'})


cfg = finputcheck(varargin,...
    {'predictAt','cell',[],{{'',[]}};
    'add_intercept','boolean',[],0;
    'plotParam','',[],{};
    'n_topos','integer',[],15,
    'baseline','real',[min(ufresult.times) max(ufresult.times)],[];...
    'betaSetName','string',fieldnames(ufresult),'beta';
    'channel','integer',[],[];
    'caxis','',[],'same';
    'figure','boolean',[0 1],1;
    },'mode','ignore');
if(ischar(cfg)); error(cfg);end

if cfg.add_intercept
    error('not supported yet')
    
end


data = ufresult.(cfg.betaSetName);

if isempty(cfg.channel)
   cfg.channel = 1:size(data,1);
end
if isempty(cfg.plotParam)
    param = 1:length(ufresult.param);
else
    param = find(ismember({ufresult.param(:).name},cfg.plotParam));
end

data = data(cfg.channel,:,param);

% Baseline Correction
if ~isempty(cfg.baseline)
   data= bsxfun(@minus,data ,mean(data(:,(ufresult.times>=cfg.baseline(1))& (ufresult.times<cfg.baseline(2)),:),2));
end

% caculate common coloraxis
if ischar(cfg.caxis) && strcmp(cfg.caxis,'same')
    cfg.caxis = prctile(data(:),[1 99]);
    cfg.caxis = [-max(abs(cfg.caxis)) max(abs(cfg.caxis))];
end

% plot the topoplots
cfg.butterfly = 'no';
ax = plot_topobutter(data,ufresult.times,ufresult.chanlocs(cfg.channel),cfg);

% plot row-names
for k = 1:length(param)
    ax.topo.topo{k}.image{1}.Units = 'normalized';
    % get the values for the rowname
    val = round(ufresult.param(param(k)).value,2,'significant');
    
    if isnan(val)
       val = []; 
    end
    str = sprintf('[%s]: %s: %g',strjoin_custom(ufresult.param(param(k)).event),ufresult.param(param(k)).name,val);
    t = text(ax.topo.topo{k}.image{1},0,1.1,10,str,'Units','normalized','HorizontalAlignment','left');
    
    % we want them on top, else they can be behind the topos
    uistack(t,'top')
    t.Interpreter = 'none'; % no interpretation
end


%% plot time-axis at the bottom
%% Define Plotting time and get an idea of what plot_topo is doing
n_topo = size(ax.topo.topo{1}.image,2);



if isfield(cfg,'time')
    timelim = cfg.time;
else
    timelim = [min(ufresult.times),max(ufresult.times)];
end

topotimes = linspace(timelim(1),timelim(2),n_topo+1);
topotimes = round(topotimes,2,'significant');

%% Define timeline axes
% We want a bit of space to the topoplots
parentPos = get(ax.topo.parent,'Position');


pos_width = parentPos(3)/(n_topo+1);

plot_pos = [parentPos(1)+ pos_width,...
    parentPos(2), ...
    parentPos(3) - pos_width*2,...line
    0];



plotAxes = axes('Position',plot_pos);
ax.timeline= plotAxes;


xlim(timelim);
set(ax.timeline,'XTick',topotimes,'YTick',[])

%%


if nargout == 1
    varargout{1} = ax;
end

