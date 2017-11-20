function dc_plotDesignmat(EEG,varargin)
%Plots the designmatrix
%If the matrix is very large (the timeexpanded/Xdc matrix) we do not plot
%everything, but only the middle 1000s. We also try to zoom in
%automatically, but this fails sometimes
%
%Arguments:
%   cfg.timeexpand' (boolean):
%        0: Plots EEG.deconv.X (default)
%        1: Plots EEG.deconv.Xdc
%   cfg.logColor(boolean): plot the color on logscale (default 0)
%   cfg.sort(boolean): Sort the designmatrix
%   cfg.figure (1/0): Open a new figure (default 1)
%
%*Example:*
% dc_plot_designmat(EEG)
%
% dc_plot_designmat(EEG,'timeexpand',1) %plot the timeexpanded X

% Secret option: 'addContData'
cfg = finputcheck(varargin,...
    {'timeexpand',   'boolean',[],0;...
    'logColor','boolean',[0,1],0;...
    'figure','boolean',[],1;...
    'sort','boolean',[0,1],0;...
    'addContData','boolean',[0,1],0;... %undocumented, adds y-data as a subplot
    },'mode','ignore');


if cfg.timeexpand
    yAxisLabel = 'time [s]';
    
    % center the plotting window around some event in middle of EEG.event
    % If we center instead around some fixed latency (e.g. middle of recording)
    % there may not be any experimental events around to see.
    % [note: this does not guarantee that the "midEvent" is one that is 
    % modeled in the linear model]
    ix_midEvent = round(length(EEG.event)/2); % take "center" event
    midEventSmp = EEG.event(ix_midEvent).latency;
    midEventLat = EEG.times(midEventSmp); % in ms
    
    % time_lim = EEG.times(ceil(end/2)) + [-100,100] * 1000;
    
    TOTALPLOTWIN_SEC = 60; % total width of default (maximum) plotting window (in sec.)
    time_lim    = midEventLat + round([-TOTALPLOTWIN_SEC/2,TOTALPLOTWIN_SEC/2]);
    
    if min(EEG.times) > time_lim(1) || max(EEG.times) < time_lim(2)
        warning('the design-matrix is too large to display, we show only the middle 60 seconds.') % THIS MSG SEEMS WRONG
    end
    
    time_ix = find(EEG.times > time_lim(1) & EEG.times < time_lim(2));
    yAxis = EEG.times(time_ix);
    X = EEG.deconv.Xdc(time_ix,:);
    shiftByOne = 0; % dont shift the XTicks by one

else
    yAxisLabel = 'event number';
    X = EEG.deconv.X;
    yAxis = 1:size(X,1);
    shiftByOne = 1;  % shift the XTicks by one
    
end
if cfg.logColor
    X = log(X);
end

if cfg.sort
    X = sortrows(X);
end

nPred = size(X,2);
nPredTheory = length(EEG.deconv.colnames);

if cfg.figure
    figure
end
if cfg.timeexpand && cfg.addContData
    subplot(1,20,[2:20])
end
ig = imagesc(1:nPred,yAxis,X);
ax = get(ig,'parent');

r = linspace(0,nPred,nPredTheory*2+1);
set(ax,'XTick',r(2+shiftByOne:2:end))
set(ax,'XTickLabel',EEG.deconv.colnames)
set(ax,'TickLabelInterpreter','none')
%set(ax,'YDir','normal')
set(ax,'YDir','reverse') % event number or time flow from top to bottom in plotted design matrix
ylabel(ax,yAxisLabel)

colorbar

if length(r)>3
    vline(r(3+shiftByOne:2:end-1),'-r')
end


if cfg.timeexpand && cfg.addContData
    subplot(1,20,1)
    plot(EEG.data(1,:),EEG.times/1000)
    %     zoom yon
    %     zoom(length(yAxis)./(size(EEG.deconv.timebasis,1)*25))
    %     pan yon
    set(gca,'YTickLabel','','XTickLabel','')
    linkaxes([ax,gca],'y' )
    
end

if cfg.timeexpand
    
    % warning('auto-zoom to useful resolution (25x the stimulus-window)')
    % warning('XXX could be dimension 2 for splines?')
    % zoom yon
    % zoom(length(yAxis)./(size(EEG.deconv.timebasis,1)*25))
    
    axes(ax)
    pan yon
    if isfield(EEG,'event')
        allTypes = {EEG.event(:).type};
        lat = [EEG.event(:).latency]/EEG.srate; % in s
        
        % look for event latencies of events in the plotted region
        latix = (lat < (time_lim(2))) & (lat > (time_lim(1)));
                    
        % plot all event onsets as horizontal lines with different colors per event
        if sum(latix) > 0 && sum(latix) < 1000 % only do this if less than 1000 events (otherwise buggy)
            
            [un,~,c] = unique(allTypes(latix));
            colorList = repmat({{'r','g','b','c','m','y','k','w'}},ceil(length(un)/8),1); % ugly and lazy... but who has more than that different types...
            colorList = [colorList{:}];
            hline(lat(latix)/1000,colorList(c))
        end
        
        % change ylim to zoom in to the event that is in middle of cut data
        %latTmp = lat(latix);
        %ylim([latTmp(ceil(end/2))/1000-3 latTmp(ceil(end/2))/1000+3])
        ZOOMPLOTWIN_SEC = 3;
        ylim([midEventLat/1000-ZOOMPLOTWIN_SEC midEventLat/1000+ZOOMPLOTWIN_SEC])
    end
    
end
