function [betaSetName] = dc_unfoldbetaSetname(unfold,varargin)
% Find out whether we want beta_dc, beta_nodc and if there are other fields
% that have the same size that we should plot as columns.

% parse inputs
cfg = finputcheck(varargin,...
    {'deconv','integer',[0,1,-1],-1, ... % -1 is autodetect
    'dataField','string','',''...
    },'mode','ignore');

nBetaSets = 1;
betaSetName = [];
if cfg.deconv == -1
    assert(isfield(unfold,'beta')|isfield(unfold,'beta_nodc'),'error: to use autodetect at least the field unfold.beta  or unfold.beta_nodc needs to exist')
    fn = fieldnames(unfold);
    
    if isfield(unfold,'beta')
        sizeBeta = size(unfold.beta);
    else
        sizeBeta = size(unfold.beta_nodc);
    end
    for f = fn'
        if strcmp(f,'times')
            continue
        end
        if length(sizeBeta) == length(size(unfold.(f{1}))) &&  all(sizeBeta == size(unfold.(f{1})))
            nBetaSets = nBetaSets+1;
            betaSetName = [betaSetName f(1)];
        end
    end
elseif cfg.deconv==0
    betaSetName = {'beta_nodc'};
else
    betaSetName = {'beta'};
    
end
