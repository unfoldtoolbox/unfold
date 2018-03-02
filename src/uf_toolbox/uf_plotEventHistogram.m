function uf_plotEventHistogram(EEG,varargin)
% Function that plots histogram of all events in the EEG.event structure
%
%Arguments:
%   cfg.eventtypes: Restrict the histogram to a specific eventtypes
%
%Return:
%
%*Example:*
% uf_plotEventHistogram(EEG,'eventA')

cfg = finputcheck(varargin,...
    {'eventtypes',   'cell', [], {};...
    },'mode','ignore');
if(ischar(cfg)); error(cfg);end

if isempty(cfg.eventtypes)
    t2 = struct2table(EEG.event);
else
    t2 = struct2table(EEG.event(ismember({EEG.event(:).type},cfg.eventtypes)));
end
escapeString = @(tStr)regexprep(tStr,'(_)','\\$1');

VarNames = escapeString(t2.Properties.VariableNames);
removeList = {'channel','bvtime','bvmknum','code','urevent'};
for k = length(VarNames):-1:1
    if any(strcmp(VarNames{k},removeList))
        VarNames(k) = [];
    end
end
figure,


for k = 1:length(VarNames)
    subplot(1,length(VarNames),k)
    data = t2{:,strcmp(escapeString(t2.Properties.VariableNames),VarNames{k})};
    if iscell(data)
        data =data(cellfun(@(x)~isempty(x),data));
        if all(cellfun(@(x)isnumeric(x),data))
            data = cell2mat(data);
        end
    end
    
    if ~isnumeric(data)
        data = categorical(data);
        [xcounts,centers] = hist(data);
        centersBar = 1:length(unique(centers));
        bar(centersBar,xcounts,'EdgeColor','none','FaceColor',[0 .5 .5])
        set(gca,'XTickLabel',centers)
    else
        
        [xcounts,centersBar] = hist(data,20); % when numeric, we cann add how many bins
        [f,xi] = ksdensity(data,[]);
        bar(centersBar,xcounts,'EdgeColor','none','FaceColor',[0 .5 .5])
        xbefore = xlim;
        hold all
        f = f*max(xcounts)/max(f);
        plot(xi,f,'Color',[0.5 0 0.5],'LineWidth',1)
        xlim(xbefore);
        
    end
    
    
    title(VarNames{k})
end

