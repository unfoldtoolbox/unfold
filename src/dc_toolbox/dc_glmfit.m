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
%   EEG set, need to have EEG.deconv.dcX compatible with
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
    { 'save_memory_or_time', 'string',{'memory','time','deprecated'}, 'deprecated';
    'method', 'string',{'par-lsmr','lsmr','matlab','pinv','glmnet'}, 'lsmr';
    'glmnetalpha','real',[],1;... # used for glmnet
    'channel','integer',[],1:EEG.nbchan;
    'debug','boolean',[],0;
    },'mode','ignore');
if(ischar(cfg)); error(cfg);end


% Backwards compatibility
switch cfg.save_memory_or_time
    case 'deprecated'
    case 'memory'
        warning('save_memory_or_time is deprecated and should be replaced by ''method''')
        cfg.method = 'lsmr';
    case 'time'
        warning('save_memory_or_time is deprecated and should be replaced by ''method''')
        cfg.method = 'matlab';
end

assert(ndims(EEG.data) ==2,'EEG.data needs to be unconcatenated, you could input it as EEG.data(:,:)')
assert(size(EEG.deconv.dcX,1) == size(EEG.data,2),'Size of designmatrix (%d,%d), not compatible with EEG data(%d,%d)',size(EEG.deconv.dcX),size(EEG.data))


    X = EEG.deconv.dcX;

disp('solving the equation system');
t = tic;
beta = nan(size(EEG.deconv.dcX,2),EEG.nbchan);

if strcmp(cfg.method,'lsmr')
    
    % go tru channels
    for e = cfg.channel
        t = tic;
        fprintf('\nsolving electrode %d (of %d electrodes in total)',e,length(cfg.channel))
        
        % use iterative solver for least-squares problems (lsmr)
        [beta(:,e),ISTOP,ITN] = lsmr(EEG.deconv.dcX,sparse(double(EEG.data(e,:)')),[],10^-8,10^-8,[],200); % ISTOP = reason why algorithm has terminated, ITN = iterations
        if ISTOP == 7
            warning(['The iterative least squares did not converge for channel ',num2str(e), ' after ' num2str(ITN) ' iterations'])
        end
        fprintf('... took %.1fs',toc(t))
        %beta(:,e) =
        %lsqr(EEG.deconv.dcX,sparse(double(EEG.data(e,:)')),[],30);
        
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
    beta = nan(size(EEG.deconv.dcX,2),EEG.nbchan);
    dcX = EEG.deconv.dcX;
    data = sparse(double(EEG.data'));
    % go tru channels
    fprintf('starting parallel loop')
    parfor e = cfg.channel
        t = tic;
        
        fprintf('\nsolving electrode %d (of %d electrodes in total)',e,length(cfg.channel))
        % use iterative solver for least-squares problems (lsmr)
        [beta(:,e),ISTOP,ITN] = lsmr(dcX,data(:,e),[],10^-10,10^-10,[],200); % ISTOP = reason why algorithm has terminated, ITN = iterations
        if ISTOP == 7
            warning(['The iterative least squares did not converge for channel ',num2str(e), ' after ' num2str(ITN) ' iterations'])
        end
        fprintf('... took %.1fs',toc(t))
        
        %beta(:,e) =
        %lsqr(EEG.deconv.dcX,sparse(double(EEG.data(e,:)')),[],30);
        
    end
    
    
elseif strcmp(cfg.method,'matlab') % save time
    
    
    if cfg.debug
        spparms('spumoni',2)
    end
    
    beta(:,cfg.channel) = EEG.deconv.dcX \ sparse(double(EEG.data(cfg.channel,:)'));
    
elseif strcmp(cfg.method,'pinv')
    Xinv = pinv(full(X));
    beta = calc_beta(EEG,Xinv,cfg.channel);
    
    
elseif strcmp(cfg.method,'glmnet')
    beta = nan(size(EEG.deconv.dcX,2)+1,EEG.nbchan); %plus one, because glmnet adds a intercept
    for e = cfg.channel
        t = tic;
        fprintf('\nsolving electrode %d (of %d electrodes in total)',e,length(cfg.channel))
        %glmnet needs double precision
        fit = cvglmnet(X,(double(EEG.data(e,:)')),'gaussian',struct('alpha',cfg.glmnetalpha));
        
        %find best cv-lambda coefficients
        beta(:,e) = cvglmnetCoef(fit,'lambda_1se')';
        
        fprintf('... took %.1fs',toc(t))

    end
    beta = beta([2:end 1],:); %put the dc-intercept last
    EEG = dc_designmat_addcol(EEG,ones(1,length(EEG.deconv.dcX)),'glmnet-DC-Correction');

    
    
end
fprintf('\n LMfit finished \n')
beta = beta'; % I prefer channels X betas (easier to multiply things to)


% We need to remove customrows, as they were not timeexpanded.
eventcell = cellfun(@(x)iscell(x(1)),EEG.deconv.eventtype)*1;
eventnan = cellfun(@(x)isnan(x(1)),EEG.deconv.eventtype(~eventcell));
eventnan = find(~eventcell);


betaOut = reshape(beta(:,1:end-length(eventnan)),size(beta,1),size(EEG.deconv.dcBasis,1),sum(~ismember(EEG.deconv.col2eventtype,eventnan)));



EEG.deconv.dcBeta = betaOut;
if length(eventnan)>0
%     EEG.betaCustomrow = beta(end+1-length(eventnan):end);
    EEG.deconv.dcBetaCustomrow = beta(:,end+1-length(eventnan):end);
end
EEG.deconv.channel = cfg.channel;

end

function [beta] = calc_beta(EEG,Xinv,channel)
beta = nan(size(Xinv,1),EEG.nbchan);
for c = channel
    beta(:,c)= (Xinv*squeeze(EEG.data(c,:,:))');
end
end