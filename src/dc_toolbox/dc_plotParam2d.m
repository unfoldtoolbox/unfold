function [] = dc_plotParam2d(t,varargin)
%Function not yet ready (sorry)
% This function plots an imagesc plot of time vs. parameter of choice
%
%Arguments:
% 'cfg.Parameter': name of the parameter AS IT OCCURS IN THE TABLE, currently
%              no search is supported to find the spline number directly.
% 'cfg.add_intercept': add the intercept to the plot, default 1
% 'cfg.type':  subset of the type variable


warning('Deprecated. This function is currently not used, but might be used to define erpimages')
cfg= finputcheck(varargin,...
    {'type','string',unique(t.type),[];
    'parameter','string',unique(t.parameter),[];
    'figure','boolean',[0 1],1;
    'add_intercept','boolean',[0 1],1; % we think interpreting results is easiest with the intercept included
    },'mode','ignore');
if(ischar(cfg))
    error(cfg);
end


select = strcmp(t.type,cfg.type);
tmpIntercept = t.signal(select&strcmp(t.parameter,'(Intercept)'));


a = t(select&strcmp(t.parameter,cfg.parameter),:);

tidx = sort(unique(a.timevec));
paramidx = unique(a.paramvalue);
b = full(sparse(arrayfun(@(x)find(x==tidx),a.timevec),arrayfun(@(x)find(x==paramidx),a.paramvalue),a.signal));

if cfg.add_intercept
    b = bsxfun(@minus,b,tmpIntercept);
end
if cfg.figure
    figure;
end
imagesc(tidx,paramidx,b')

% Design Stuff
title(cfg.type)
ylabel(cfg.parameter)
xlabel('time [s]')
set(gca,'YDir','normal')
colorbar
text(1,-0.1,['Intercept Added: ' num2str(cfg.add_intercept)],'units','normalized')