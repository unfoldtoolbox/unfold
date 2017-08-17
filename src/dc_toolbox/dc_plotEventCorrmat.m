function corrOut = dc_plotEventCorrmat(EEG,varargin)
% Plots a correlation matrix of the event structure
%
%Arguments:
%   eventtype (cell): Subselect the eventtype, by default chooses all
%   figure (0/1): whether the corrmat should be plotted (default) or only
%               returned
%
%Return:
%   correlationMatrix


cfg = finputcheck(varargin,...
    {'eventtype',   'cell', [], {};...
    'figure','boolean',[],1;...
    },'mode','ignore');
if(ischar(cfg)); error(cfg);end

if isempty(cfg.eventtype)
    t2 = struct2table(EEG.event);
else
    t2 = struct2table(EEG.event(ismember({EEG.event(:).type},cfg.eventtype)));
end
if size(t2,1) == 0 && isempty(cfg.eventtype)
    error('EEG.event seems to be empty')
elseif   size(t2,1) == 0 && ~isempty(cfg.eventtype)
    error('EEG.event is empty for event(s): %s',cfg.eventtype{:})
end

t2 = t2(:,varfun(@isnumeric,t2,'output','uniform'));
t2.time = t2.latency/EEG.srate;
t2.fixDurCurrent = diff([t2.time;nan(1)]);
t2.fixDurPrevious = diff([nan(1);t2.time]);

nanrows = any(isnan(t2{:,:}),2);
corrData = corr(t2{~nanrows,:});
if cfg.figure
    figure,
    imagesc(corrData)
    set(gca,'YTick',1:length(corrData))
    set(gca,'YTickLabel',t2.Properties.VariableNames)
    set(gca,'XTick',1:length(corrData))
    set(gca,'XTickLabel',t2.Properties.VariableNames)
    
    
    caxis([-1 1]),colorbar
end
corrOut = array2table(corr(t2{~isnan(t2.fixDurCurrent) & ~isnan(t2.fixDurPrevious),:}),'VariableNames',t2.Properties.VariableNames,'RowNames',t2.Properties.VariableNames);