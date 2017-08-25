function [] = dc_plotParam2d(unfold,varargin)
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
    {'plotParam','cell',[],{};
    'add_intercept','boolean',[0 1],0;
    'channel','integer',[],[];
    },'mode','ignore');
if(ischar(cfg))
    error(cfg);
end

if cfg.add_intercept
    error('adding intercept not yet supported')
end

paramIdx = find(strcmp(unfold.deconv.variableType,'spline') | strcmp(unfold.deconv.variableType,'continuous'));
if cfg.add_intercept
    b = bsxfun(@minus,b,tmpIntercept);
end
betaSetName = 'beta_nodc'
if any(strcmp(unfold.deconv.variableType,'spline')) && size(unfold.(betaSetName),3) > size(unfold.deconv.predictorSplines{1}.spline2val,2)

end

for p = 1:length(paramIdx)
    figure
    varName = unfold.deconv.variableNames(paramIdx(p));
    unfoldTmp = dc_getParam(unfold);
    ix = strcmp(varName,{unfold.epoch.name});
    val = [unfold.epoch.value];
    imagesc(unfold.times,val(ix),squeeze(unfold.(betaSetName)(cfg.channel,:,ix))')
    
    ylabel(varName)
    xlabel('time [s]')
    set(gca,'YDir','normal')
    colorbar
    text(1,-0.1,['Intercept Added: ' num2str(cfg.add_intercept)],'units','normalized')
    
end

% Design Stuff

