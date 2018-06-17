function [] = uf_plotParam2d(ufresult,varargin)
% Plots a 2D plot of parameter vs. time with the predicted value as the
% third dimension
% This function plots an imagesc plot of time vs. parameter of choice
%
%Arguments:
% 'plotParam'     : Name of parameter to be plotted. can be empty to plot
%                   all splines/continuous parameters
% 'add_intercept' : add the intercept to the plot, default 0
% 'channel'       : Specify which channel-idx to plot
% 'betaSetName'   : Default 'beta'. Can be any field of the ufresult-struct
% 'caxis'         : Default [], specify coloraxis
%
%Example:
%uf_plotParam2d('plotParam','continuosPredictorA')

cfg= finputcheck(varargin,...
    {'plotParam','',[],{};
    'betaSetName','string',fieldnames(ufresult),'beta';
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

% check whether the user tried to enter EEG.unfold directly into this 
% function without running uf_condense first
if ~isfield(ufresult,'param') & isfield(ufresult,'unfold')
    error('\n%s(): You cannot directly enter the unfold output into this function - you have to run uf_condense() first',mfilename)
end

paramIdx = find(strcmp(ufresult.unfold.variabletypes,'spline') | strcmp(ufresult.unfold.variabletypes,'continuous'));

for p = 1:length(paramIdx)
    
    if any(strcmp(ufresult.unfold.variabletypes,'continuous'))
        cfgtmp= struct();
        cfgtmp.auto_n = 100;
        cfgtmp.auto_method = 'linear';
        ufresult = uf_predictContinuous(ufresult,cfgtmp);
    end
    figure
    varName = ufresult.unfold.variablenames(paramIdx(p));
    if ~isempty(cfg.plotParam) && sum(strcmp(varName,cfg.plotParam))==0
        continue
    end
    ix = strcmp(varName,{ufresult.param.name});
    val = [ufresult.param.value];
    
    b = squeeze(ufresult.(cfg.betaSetName)(cfg.channel,:,:));
    if cfg.add_intercept
        % get eventtname of current predictor
        evtType =  ufresult.unfold.cols2eventtypes(paramIdx);
        % combine it in case we have multiple events
        evtType = strjoin_custom(ufresult.unfold.eventtypes{evtType},'+');
        % find the events in ufresult.param
        epochStr = cellfun(@(x)strjoin_custom(x,'+'),{ufresult.param(:).event},'UniformOutput',0);
        % check the event to be the same & that it is an intercept
        interceptIdx = strcmp(evtType,epochStr) & strcmp('(Intercept)',{ufresult.param(:).name});
        
        if sum(interceptIdx) ~= 1
            error('could not find the intercept you requested (predictor:%s',varName)
        end
        b(:,ix) = bsxfun(@plus,b(:,ix),b(:,interceptIdx));
    end
    
    
    imagesc(ufresult.times,val(ix),b(:,ix)')
    
    ylabel(varName,'Interpreter','none')
    xlabel('time [s]')
    if isfield(ufresult,'chanlocs')&&~isempty(ufresult.chanlocs)
        chanName = ufresult.chanlocs(1).labels;
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
vline(0,':k')
    
    
end
