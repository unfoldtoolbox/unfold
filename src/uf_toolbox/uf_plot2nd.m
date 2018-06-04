function ax = uf_plot2nd(d2nd,varargin)
%Preliminary function to plot multiple subjects (2nd-level analysis)
% This function allows to plot multiple subjects at the same time
% the function requires the data to be in the following format:
% ufresult.beta(CHAN,TIME,PARAM,SUBJECT)
% 
% Each line is one subject, its possible to calculate confidence intervals
%
%Arguments:
%   cfg.channel: Which channel to use
%       
%   cfg.plotParam: (default 1, as in glmnet), can be 0 for L2 norm, 1 for L1-norm or
%       something inbetween for elastic net
%
%   cfg.bootci:  (default 1) calculate and plot boostraped confidence intervals
%
%   cfg.singlesubjects: (default 1) plot the singlesubject lines
%
%   ... :      Other parameters are linked to uf_plotParam
%   
%Returns:
%   nothing
%
%Example:
%   uf_plot2nd(ufresult2nd,'channel',2)
%
cfg = finputcheck(varargin,...
    {'channel','integer',[],[];
    
    'plotParam','',[],[];
    'bootci','boolean',[],1;
    'singlesubjects','boolean',[],1;
    
    },[],'ignore');
if(ischar(cfg)); error(cfg);end

assert(isstruct(d2nd),'input need to be a struct of subjects')

% first needs a check whether beta or beta_nodc is present
% assert(size(beta,4) ~= 1,'multiple subjects need to be present')

% TODO: preselct by PlotParam
if any(cellfun(@(x)~isempty(x),strfind({d2nd.param.type},'spline')))
    cfg.withSpline = 1;
else
    cfg.withSpline = 0; %allows for faster method
end

if cfg.withSpline && cfg.bootci
    warning('currently cannot plot bootci with converted splines')
end

cfgPlot = cfg;
% cfgPlot.add_intercept = 0;
if isempty(cfg.plotParam)
    cfgPlot = rmfield(cfgPlot,'plotParam');
else
    cfgPlot.plotParam = cfg.plotParam;
end

if cfg.singlesubjects
    
    if ~cfg.withSpline
        d2nd2 = d2nd;
        
        d2nd2.param = repmat(d2nd2.param,1,size(d2nd.beta,4));
        % d2nd2.unfold.X = ones(1,size(d2nd2.param,2));
        if isfield(d2nd,'beta_nodc')
            d2nd2.beta_nodc = d2nd2.beta_nodc(:,:,:);
        end
        d2nd2.beta = d2nd2.beta(:,:,:);
        
        g = uf_plotParam(d2nd2,cfgPlot);
        
        
        % Make the lines black & transparent
        for a = g.facet_axes_handles(:)'
            a = get(a,'Children');
            for child = a'
                if child.LineStyle == '-'
                    child.Color = [0 0 0 0.1];
                end
            end
        end
        
        % add mean + bootCI
        if cfg.bootci
            update(g);
            stat_summary(g,'geom','area','type','bootci');
            draw(g);
            
        end
        
    else % withspline
        %% plot unfold vs nodeconv
        % this works with splines but is slooooww
        
        for s = 1:size(d2nd.beta,4)
            d2nd2 = d2nd;
            d2nd2.beta = d2nd2.beta(:,:,:,s);
            d2nd2.beta_nodc = d2nd2.beta_nodc(:,:,:,s);
            if s >1
                cfg.gramm = g;
            end
            g = uf_plotParam(d2nd2,cfgPlot);
        end
        
        
        for a = g.facet_axes_handles(:)'
            a = get(a,'Children');
            for child = a'
                if child.LineStyle == '-'
                    child.Color = [0 0 0 0.1];
                end
            end
        end
    end
end
if cfg.singlesubjects
    cfgPlot.gramm = g;
end
d2nd2 = d2nd;
d2nd2.beta = mean(d2nd2.beta(:,:,:,:),4);
d2nd2.beta_nodc = mean(d2nd2.beta_nodc(:,:,:,:),4);
if length(d2nd2.unfold) >1
    d2nd2.unfold  = d2nd2.unfold(1);
end
ax = uf_plotParam(d2nd2,cfgPlot);

