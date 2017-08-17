function [varargout] = dc_plotParamTopo(unfold,varargin)
%Generates topoplots over time. each row is a predictor or predictor-value-combination
%
%Arguments:
% 'plotParam' : cell array of parameters to be plotted, if empty plots all
% 'pred_value': only used for continuous / spline predictors. 
%               a cell of cell arrays, e.g. {{'parName',linspace(0,10,5)},{'parname2',1:5}}
%               This would split up the parName-predictor into 5 bins from
%               0 to 10, so 5 rows of topoplots would be plotted. Default are 7 lines
%               from min to max (see code for more info)
%
% 'add_intercept' : (not yet supported) whether to add the intercept to each curve (thus not
%               the 'pure' effect of a independent variable is plotted but
%               more of an ERP-like plot is generated
%
% 'chan' : plot only a subset of channels
%
% 'caxis' ('same',default:[]) if 'same', generates the same coloraxis based
%       on the 95% percentile of the selected beta-values. can be
%       customized to whichever caxis
% 'figure'    : plot in new figure (1) or old (0), default: (1)
%
%Returns:
%   structure of all plotting axes.
%
%*Examples:*
%  dc_plotParamTopo(EEG,'plotParam',{'FactorX','FactorC'})


cfg = finputcheck(varargin,...
    {'pred_value','cell',[],{{'',[]}};
    'add_intercept','boolean',[],0;
    'plotParam','cell',[],{};
    'channel','integer',[],[];
    'caxis','',[],'same';
    'figure','boolean',[0 1],1;
    },'mode','ignore');
if(ischar(cfg)); error(cfg);end

% data_bsl = bsxfun(@minus,unfold.data,mean(unfold.data(:,unfold.times<-0.1,:),2));
% plot_topobutter(data_bsl(1:64,:,:),unfold.times,unfold.chanlocs(1:64))

if isempty(cfg.channel)
   cfg.channel = 1:size(unfold.beta,1);
end
if isempty(cfg.plotParam)
    param = 1:length(unfold.epoch);
else
    param = find(ismember(cfg.plotParam,{unfold.epoch(:).name}));
end

unfold.beta = unfold.beta(cfg.channel,:,param);
if ischar(cfg.caxis) && strcmp(cfg.caxis,'same')
    cfg.caxis = prctile(unfold.beta(:),[1 99]);
    cfg.caxis = [-max(abs(cfg.caxis)) max(abs(cfg.caxis))];
end

cfg.butterfly = 'no';
ax = plot_topobutter(unfold.beta,unfold.times,unfold.chanlocs(cfg.channel),cfg);


for k = 1:length(param)
    ax.topo.topo{k}.image{1}.Units = 'normalized';
    val = round(unfold.epoch(param(k)).value,2,'significant');
    if isnan(val)
       val = []; 
    end
    t = text(ax.topo.topo{k}.image{1},0,1.1,10,sprintf('[%s]: %s: %g',strjoin(unfold.epoch(param(k)).event),unfold.epoch(param(k)).name,val),'Units','normalized','HorizontalAlignment','left');
    uistack(t,'top')
    t.Interpreter = 'none';
end

% plot time-axis at the bottom


%% Define Plotting time and get an idea of what plot_topo is doing
n_topo = size(ax.topo.topo{1}.image,2);

topotimes = linspace(min(unfold.times),max(unfold.times),n_topo+1);
round(topotimes,2,'significant')

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



xlim([min(unfold.times),max(unfold.times)]);
set(ax.timeline,'XTick',topotimes,'YTick',[])

%%


if nargout == 1
    varargout{1} = ax;
end

