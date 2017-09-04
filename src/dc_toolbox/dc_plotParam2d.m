function [] = dc_plotParam2d(unfold,varargin)
%Function not yet ready (sorry)
% This function plots an imagesc plot of time vs. parameter of choice
%
%Arguments:
% 'plotParam'     : Name of parameter to be plotted. can be empty to plot
%                   all splines/continuous parameters
% 'add_intercept' : add the intercept to the plot, default 0
% 'channel'       : Specify which channel-idx to plot
% 'betaSetName'   : Default 'beta'. Can be any field of the unfold-struct
% 'caxis'         : Default [], specify coloraxis


cfg= finputcheck(varargin,...
    {'plotParam','',[],{};
    'betaSetName','string',fieldnames(unfold),'beta';
    'add_intercept','boolean',[0 1],0;
    'caxis','',[],[];
    'channel','integer',[],[];
    },'mode','ignore');
if(ischar(cfg))
    error(cfg);
end


if isempty(cfg.channel)
    error('you need to specify a channel')
end


paramIdx = find(strcmp(unfold.deconv.variableType,'spline') | strcmp(unfold.deconv.variableType,'continuous'));

for p = 1:length(paramIdx)
    
    if any(strcmp(unfold.deconv.variableType,'continuous'))
        cfg.auto_n = 100;
        cfg.auto_method = 'linear';
        unfold = dc_getParam(unfold,cfg);
    end
    figure
    varName = unfold.deconv.variableNames(paramIdx(p));
    if ~isempty(cfg.plotParam) && sum(strcmp(varName,cfg.plotParam))==0
        continue
    end
    ix = strcmp(varName,{unfold.epoch.name});
    val = [unfold.epoch.value];
    
    b = squeeze(unfold.(cfg.betaSetName)(cfg.channel,:,:));
    if cfg.add_intercept
        % get eventtname of current predictor
        evtType =  unfold.deconv.col2eventtype(paramIdx);
        % combine it in case we have multiple events
        evtType = strjoin(unfold.deconv.eventtype{evtType},'+');
        % find the events in unfold.epoch
        epochStr = cellfun(@(x)strjoin(x,'+'),{unfold.epoch(:).event},'UniformOutput',0);
        % check the event to be the same & that it is an intercept
        interceptIdx = strcmp(evtType,epochStr) & strcmp('(Intercept)',{unfold.epoch(:).name});
        
        if sum(interceptIdx) ~= 1
            error('could not find the intercept you requested (predictor:%s',varName)
        end
        b(:,ix) = bsxfun(@plus,b(:,ix),b(:,interceptIdx));
    end
    
    
    imagesc(unfold.times,val(ix),b(:,ix)')
    
    ylabel(varName,'Interpreter','none')
    xlabel('time [s]')
    if isfield(unfold,'chanlocs')&&~isempty(unfold.chanlocs)
        chanName = unfold.chanlocs(1).labels;
    else
        chanName = num2str(cfg.channel);
    end
    set(gcf,'Name',sprintf('Unfold-Toolbox Channel: %s',chanName))
    set(gca,'YDir','normal')
    colorbar
    
    if ~isempty(cfg.caxis)
        caxis(cfg.caxis)
        
    else
        c = max(abs(caxis));
        caxis([-c c])
    end
    text(1,-0.05,sprintf('Intercept added?: %i \nBetaSetName: %s',cfg.add_intercept,cfg.betaSetName),'units','normalized')
    
    
end
