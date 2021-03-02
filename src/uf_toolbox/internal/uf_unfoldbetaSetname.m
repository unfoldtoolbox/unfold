function [betaSetName] = uf_unfoldbetaSetname(ufresult,varargin)
% Find out whether we want beta_dc, beta_nodc and whether there are other fields
% that have the same size that we should plot as columns

% parse inputs
cfg = finputcheck(varargin,...
    {'deconv','integer',[0,1,-1],-1; ... % -1 is autodetect
    'dataField','string','',''...
    },'mode','ignore');

nBetaSets     = 1;
betaSetName   = [];

if cfg.deconv == -1

    % autodetect 
    assert(~isempty(cfg.dataField) || isfield(ufresult,'beta') || isfield(ufresult,'beta_nodc'),'Error: to use autodetect at least the field ufresult.beta or ufresult.beta_nodc needs to exist')
    fn = fieldnames(ufresult);
    
    % get size of betas
    if ~isempty(cfg.dataField)
        sizeBeta = size(ufresult.(cfg.dataField));
    elseif isfield(ufresult,'beta')
        sizeBeta = size(ufresult.beta);
    else
        sizeBeta = size(ufresult.beta_nodc);
    end
    
    for f = fn'
        % check for field "times"
        if strcmp(f,'times')
            continue
        end
        
        %
        if length(sizeBeta) == length(size(ufresult.(f{1}))) &&  all(sizeBeta == size(ufresult.(f{1})))
            nBetaSets   = nBetaSets+1;
            betaSetName = [betaSetName f(1)];
        end
    end
    
elseif cfg.deconv==0
    % massive univariate model
    betaSetName = {'beta_nodc'};
else
    % deconvolution model
    betaSetName = {'beta'};
end
