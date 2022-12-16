function uf_plotDesignmat(EEG,varargin)
%% Plots the designmatrix
%If the matrix is very large (the timeexpanded/Xdc matrix) we do not plot
%everything, but only the middle 60s. In addition (for timeexpand) we plot
%the events as horizontal lines.
%
%Arguments:
%   cfg.timeexpand' (boolean):
%        0: Plots EEG.unfold.X (default)
%        1: Plots EEG.unfold.Xdc
%   cfg.logColor(boolean): plot the color on logscale (default 0)
%   cfg.sort(boolean): Sort the designmatrix (only possible for X, not Xdc)
%   cfg.figure (1/0): Open a new figure (default 1)
%   cfg.totalplotwin_sec (real): How many seconds to plot designmatrix (default: 60s)
%
%*Example:*
% uf_plot_designmat(EEG)
% uf_plot_designmat(EEG,'sort',1)
% uf_plot_designmat(EEG,'timeexpand',1) %plot the timeexpanded X

% Secret option: 'addContData'
cfg = finputcheck(varargin,...
    {'timeexpand',   'boolean',[],0;...
    'logColor','boolean',[0,1],0;...
    'figure','boolean',[],1;...
    'sort','boolean',[0,1],0;...
    'totalplotwin_sec','real',[],60;... % total width of default (maximum) plotting window (in sec.)
    'addContData','boolean',[0,1],0;... %undocumented, adds y-data as a subplot
    },'mode','error');
if ischar(cfg)
    error(cfg)
end
assert(~(cfg.timeexpand & cfg.sort),'cannot plot Xdc sorted')

if cfg.timeexpand
    yAxisLabel = 'time [s]';
    
    % center the plotting window around some event in middle of EEG.event
    % If we center instead around some fixed latency (e.g. middle of recording)
    % there may not be any experimental events around to see.
    
    
    modeledEvents = [EEG.unfold.eventtypes{:}];
    % we remove the TRFs
    eventNan = cellfun(@(x)all(isnan(x)),modeledEvents);
    modeledEvents(eventNan) = [];
    
    % temporary change "nan" to "trf" for later plotting label
    
    EEG.unfold.eventtypes(eventNan) = {{'trf'}};
    
    if isempty(modeledEvents)
        % only TRFs
        midEventLat = EEG.times(round(end/2));
        
    else
        eventstruct = EEG.event;
        eventstruct(~ismember({eventstruct.type},modeledEvents)) = [];
        ix_midEvent = round(length(eventstruct)/2); % take "center" event
        midEventSmp = round(eventstruct(ix_midEvent).latency);
        midEventLat = EEG.times(midEventSmp); % in ms
    end
    
    % time_lim = EEG.times(ceil(end/2)) + [-100,100] * 1000;
    
    time_lim    = midEventLat + round([-cfg.totalplotwin_sec/2,cfg.totalplotwin_sec/2]*1000);
    
    time_ix = find(EEG.times > time_lim(1) & EEG.times < time_lim(2));
    yAxis = EEG.times(time_ix)/1000;
    X = EEG.unfold.Xdc(time_ix,:);
    shiftByOne = 0; % dont shift the XTicks by one
    
else
    yAxisLabel = 'event number';
    X = EEG.unfold.X;
    yAxis = 1:size(X,1);
    shiftByOne = 1;  % shift the XTicks by one
    
end
if cfg.logColor
    Xpos = X(:)>0;
    Xneg = X(:)<0;
    X(Xpos) = log(X(Xpos));
    X(Xneg) = -log(-X(Xneg));
end

if cfg.sort
    X = sortrows(X);
end

nPred = size(X,2);
nPredTheory = length(EEG.unfold.colnames);

if cfg.figure
    figure
end
if cfg.timeexpand && cfg.addContData
    subplot(1,20,[2:20])
end
ig = imagesc(1:nPred,yAxis,X);
ax = get(ig,'parent');

%custom colormap
cmap = cbrewer('div','RdBu',100);
cmap=cmap(end:-1:1,:); %to have blue negative, red pos
colormap(cmap);

%force symmetric colorscale
cl = get(gca,'clim');
cl = max(abs(cl));
set(gca,'clim',[-cl cl])

r = linspace(0,nPred,nPredTheory*2+1);
set(ax,'XTick',r(2+shiftByOne:2:end))

eventlabtmp = cellfun(@(x)['evt:' strjoin(x,'+')],EEG.unfold.eventtypes,'UniformOutput',0);
eventlabels = eventlabtmp(EEG.unfold.cols2eventtypes);
predlabels = EEG.unfold.colnames;

% escape underscores
escapeString = @(tStr)regexprep(tStr,'(_)','\\$1');

% concatenate event and predictorlabel
xlab = cellfun(@(pred,evt)[escapeString(evt) '\newline' escapeString(pred)],predlabels,eventlabels,'UniformOutput',0);


set(ax,'XTickLabel',xlab)
set(ax,'XTickLabelRotation',45);
set(ax,'YDir','reverse') % event number or time flow from top to bottom in plotted design matrix
set(ax,'box','off')

ylabel(ax,yAxisLabel)

colorbar

if length(r)>3
    vline(r(3+shiftByOne:2:end-1),'-k')
end


if cfg.timeexpand && cfg.addContData
    subplot(1,20,1)
    plot(EEG.data(1,:),EEG.times/1000)
    %     zoom yon
    %     zoom(length(yAxis)./(size(EEG.unfold.timebasis,1)*25))
    %     pan yon
    set(gca,'YTickLabel','','XTickLabel','')
    linkaxes([ax,gca],'y' )
    
end

if cfg.timeexpand
    
    % warning('auto-zoom to useful resolution (25x the stimulus-window)')
    % warning('XXX could be dimension 2 for splines?')
    % zoom yon
    % zoom(length(yAxis)./(size(EEG.unfold.timebasis,1)*25))
    
    axes(ax)
    pan yon
    if isfield(EEG,'event')
        allTypes = {EEG.event(:).type};
        lat = [EEG.event(:).latency]/EEG.srate*1000; % in ms
        
        % look for event latencies of events in the plotted region
        latix = find((lat < (time_lim(2))) & (lat > (time_lim(1))));
        
        % plot all event onsets as horizontal lines with different colors per event
        if ~isempty(latix) && length(latix) < 1000 % only do this if less than 1000 events (otherwise buggy)
            
            [un,~,c] = unique(allTypes(latix));
            makeGray = ~ismember(un,modeledEvents);
            %colorList = cbrewer('qual','Set1',max(length(un),3));
            colorList = cbrewer('qual','Dark2',max(length(un),3));
            legendlist = {};
            for ev = 1:length(un)
                eventlat = lat(latix(c==ev))/1000;
                if isempty(eventlat)
                    continue
                end
                if ismember(ev,find(makeGray))
                    col = [0.8 0.8 0.8];
                else
                    col = colorList(ev,:);
                end
                lineProp = {'color',col,'linewidth',0.7}; % otherwise lines hard to see
                if length(eventlat) > 1
                    lineProp = repmat({lineProp},length(eventlat),1);
                end
                
                hh = hline(eventlat,'','',lineProp);
                legendlist{end+1,1} = un{ev};
                legendlist{end,2} = ismember(ev,find(makeGray));
                legendlist{end,3} = hh(1);
                
                
            end
            if sum(makeGray)>0
                legendlist{end+1,1} = 'non-modeled';
                legendlist{end,2} = 0;
                legendlist{end,3} = legendlist{find([legendlist{:,2}],1),3};
                legendlist(find([legendlist{:,2}]),:) = [];
                
                
            end
            legend([legendlist{:,3}],[legendlist(:,1)])
        else
            warning('more than 1000 events found in the restricted timeframe, this makes Matlab crash and we are not plotting them')
        end
        
                
        % change ylim to zoom in to the event that is in middle of cut data
        %latTmp = lat(latix);
        %ylim([latTmp(ceil(end/2))/1000-3 latTmp(ceil(end/2))/1000+3])
        ZOOMPLOTWIN_SEC = 3;
        ylim([midEventLat/1000-ZOOMPLOTWIN_SEC midEventLat/1000+ZOOMPLOTWIN_SEC])
        fprintf('Zoomed in version (~%is) of time-expanded designmatrix Xdc (shown are the middle %is). Use ''hand''-tool to drag the plot up/down. \n',2*ZOOMPLOTWIN_SEC,cfg.totalplotwin_sec)
    end
    
end
