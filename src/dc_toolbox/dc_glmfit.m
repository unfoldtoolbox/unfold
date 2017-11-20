function [EEG,beta] = dc_glmfit(EEG,varargin)
%% Fit the fullX designmatrix on the data and returns beta and stats
% This function solves the Equation X*beta = EEG.data, with X = Designmat.
% There are multiple algorithms implemented, a slow iterative algorithm
% that runs on sparse matrices (default) that solves each channel in turn
% and the matlab algorithm which solves all channels at the same time, but
% take quite a lot of memory.
%
%Arguments:
%   cfg.method (string):
%    * "lsmr"      default; SAVES MEMORY an iterative solver is used, this is
%    very memory efficient, but is a lot slower than the 'time' option
%    because each electrode has to be solved independently. The LSMR
%    algorithm is used for sparse iterative solving.
%   
%    * "par-lsmr"  same as lsmr, but uses parfor with ncpu-1.
%
%    * "matlab"    , uses matlabs native A/b solver. For moderate to big
%    design-matrices it will need *a lot* of memory (40-60GB is easily
%    reached)
%
%    * "pinv"      A naive pseudo-inverse, generally not recommended due to
%    floating point instability
%
%    * "glmnet"    uses glmnet to fit the linear system. This by default uses
%    L1-Norm aka lasso (specified as cfg.glmnetalpha = 1). For
%    ridge-regression (L2-Norm) use (cfg.glmnetalpha = 0). Something
%    inbetween results in elastic-net. We use the cvglmnet functionality
%    that automatically does crossvalidation to estimate the lambda
%    parameter (i.e. how strongly parameter values should be regularised
%    compared to the fit of the model). We use the glmnet recommended
%    'lambda_1se', i.e. minimum lambda + 1SE buffer towards more strict
%    regularisation.
%
%   cfg.lsmriterations: (default 200), defines how many steps the iterative
%                   solver should search for a solution. 
%
%   cfg.glmnetalpha: (default 1), can be 0 for L2 norm, 1 for L1-norm or
%                    something inbetween for elastic net
%
%   cfg.channel(array): Restrict the beta-calculation to a subset of
%               channels. Default is all channels
%
%   cfg.debug (boolean): 0, only with method:matlab, outputs additional
%                  details from the solver used
%
%   save_memory_or_time: Deprecated, has been renamed to 'method' EEG:  the
%   EEG set, need to have EEG.deconv.Xdc compatible with
%         the size of EEG.data
%
%Return:
% EEG.deconv.beta: array (nchan x ntime x npred) (ntime could be
% n-timesplines, n-fourierbasis or samples)
%
%*Examples:*
% EEG = dc_glmfit(EEG); EEG = dc_glmfit(EEG,'method','matlab');
%

fprintf('\ndc_glmfit(): Fitting deconvolution model...');


cfg = finputcheck(varargin,...
    {'method', 'string',{'par-lsmr','lsmr','matlab','pinv','glmnet'}, 'lsmr';
    'lsmriterations','integer',[],200;
    'glmnetalpha','real',[],1;... # used for glmnet
    'channel','integer',[],1:size(EEG.data,1);
    'debug','boolean',[],0;
    },'mode','ignore');
if(ischar(cfg)); error(cfg);end




assert(ndims(EEG.data) ==2,'EEG.data needs to be unconcatenated. Did you epoch your data already? We need continuous data for this fit')
assert(size(EEG.deconv.Xdc,1) == size(EEG.data,2),'Size of designmatrix (%d,%d), not compatible with EEG data(%d,%d)',size(EEG.deconv.Xdc),size(EEG.data))
assert(~any(isnan(EEG.deconv.Xdc(:))),'Warning NAN values found in designmatrix. will not continue')

X = EEG.deconv.Xdc;

disp('solving the equation system');
t = tic;
beta = nan(size(EEG.deconv.Xdc,2),EEG.nbchan);
data = EEG.data;

%% Remove data that is unnecessary for the fit
% this helps calculating better tolerances for lsmr
emptyRows = sum(abs(X),2) == 0;
X(emptyRows,:)  = [];
data(:,emptyRows) = [];

% % normalize data
% READ THIS
% this is uncommented because it doesn't actually solve the problem I
% wanted to solve. nevertheless it might be helpful to normalize data
% before the fits. We hav eto check the methods, some do this automatically
% (glmnet). Others might (matlab?). Pinv should benefit from this, lsmr - I
% don't really know.
%
% if any(strcmp(cfg.method,{'lsmr','par-lsmr'}))
%     % I don't do mean/median subtraction yet, because I'm not actually sure
%     % how to undo it again? I don't have a total-EEG signal Intercept.
%     % Maybe I could just do it and not undo it again?
%
%     normalize = struct();
%     %     normalize.median = nanmedian(data,2);
%     normalize.mad = mad(data')';
%     %     data = bsxfun(@minus,data,normalize.median);
%     data = bsxfun(@rdivide,data,normalize.mad);
% end



if strcmp(cfg.method,'lsmr')
    
    % go tru channels
    for e = cfg.channel
        t = tic;
        fprintf('\nsolving electrode %d (of %d electrodes in total)',e,length(cfg.channel))
        
        % use iterative solver for least-squares problems (lsmr)
        [beta(:,e),ISTOP,ITN] = lsmr(X,sparse(double(data(e,:)')),[],10^-8,10^-8,[],cfg.lsmriterations); % ISTOP = reason why algorithm has terminated, ITN = iterations
        if ISTOP == 7
            warning(['The iterative least squares did not converge for channel ',num2str(e), ' after ' num2str(ITN) ' iterations'])
        end
        fprintf('... %i iterations, took %.1fs',ITN,toc(t))
        %beta(:,e) =
        %lsqr(EEG.deconv.Xdc,sparse(double(EEG.data(e,:)')),[],30);
        
    end
    
elseif strcmp(cfg.method,'par-lsmr')
    fprintf('starting parpool with ncpu-1...')
    pools = gcp('nocreate');
    cpus = feature('numCores');
    if size(pools) == 0
        parpool(cpus-1);
    end
    fprintf('done\n')
    addpath('../lib/lsmr/')
    beta = nan(size(X,2),EEG.nbchan);
    Xdc = X;
    data = sparse(double(data'));
    % go tru channels
    fprintf('starting parallel loop')
    parfor e = cfg.channel
        t = tic;
        
        fprintf('\nsolving electrode %d (of %d electrodes in total)',e,length(cfg.channel))
        % use iterative solver for least-squares problems (lsmr)
        [beta(:,e),ISTOP,ITN] = lsmr(Xdc,data(:,e),[],10^-10,10^-10,[],200); % ISTOP = reason why algorithm has terminated, ITN = iterations
        if ISTOP == 7
            warning(['The iterative least squares did not converge for channel ',num2str(e), ' after ' num2str(ITN) ' iterations'])
        end
        fprintf('... took %.1fs',toc(t))
        
        %beta(:,e) =
        %lsqr(EEG.deconv.Xdc,sparse(double(EEG.data(e,:)')),[],30);
        
    end
    
    
elseif strcmp(cfg.method,'matlab') % save time
    
    
    if cfg.debug
        spparms('spumoni',2)
    end
    
    beta(:,cfg.channel) = X \ sparse(double(data(cfg.channel,:)'));
    
elseif strcmp(cfg.method,'pinv')
    Xinv = pinv(full(X));
    beta = calc_beta(data,EEG.nbchan,Xinv,cfg.channel);
    
    
elseif strcmp(cfg.method,'glmnet')
    beta = nan(size(X,2)+1,EEG.nbchan); %plus one, because glmnet adds a intercept
    for e = cfg.channel
        t = tic;
        fprintf('\nsolving electrode %d (of %d electrodes in total)',e,length(cfg.channel))
        %glmnet needs double precision
        fit = cvglmnet(X,(double(data(e,:)')),'gaussian',struct('alpha',cfg.glmnetalpha));
        
        %find best cv-lambda coefficients
        beta(:,e) = cvglmnetCoef(fit,'lambda_1se')';
        fit.glmnet_fit = [];
        EEG.deconv.glmnet(e) = fit;
        fprintf('... took %.1fs',toc(t))
        
    end
    beta = beta([2:end 1],:); %put the dc-intercept last
    EEG = dc_designmat_addcol(EEG,ones(1,size(EEG.deconv.Xdc,1)),'glmnet-DC-Correction');
    
    
    
end
fprintf('\n LMfit finished \n')
beta = beta'; % I prefer channels X betas (easier to multiply things to)

% Unnormalize data
% if any(strcmp(cfg.method,{'lsmr','par-lsmr'}))
%    beta = bsxfun(@times,beta,normalize.mad);
%
% end

% We need to remove customrows, as they were not timeexpanded.
eventcell = cellfun(@(x)iscell(x(1)),EEG.deconv.eventtype)*1;
eventnan = cellfun(@(x)isnan(x(1)),EEG.deconv.eventtype(~eventcell));
eventnan = find(~eventcell);


betaOut = reshape(beta(:,1:end-length(eventnan)),size(beta,1),size(EEG.deconv.timebasis,1),sum(~ismember(EEG.deconv.col2eventtype,eventnan)));



EEG.deconv.beta_dc = betaOut;
if length(eventnan)>0
    %     EEG.betaCustomrow = beta(end+1-length(eventnan):end);
    EEG.deconv.beta_dcCustomrow = beta(:,end+1-length(eventnan):end);
end
EEG.deconv.channel = cfg.channel;

end

function [beta] = calc_beta(data,nbchan,Xinv,channel)
beta = nan(size(Xinv,1),nbchan);
for c = channel
    beta(:,c)= (Xinv*squeeze(data(c,:,:))');
end
end