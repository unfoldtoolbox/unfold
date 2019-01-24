function [EEG] = uf_glmfit_nodc(EEG,varargin)
%% A function to solve the inverse problem without deconvolution
% Simple function to do massive univariate linear model. The function
% expects EEG.data to be (CHAN,TIME,EPOCH) with EPOCH the same number as
% EEG.unfold.X. 
%
% It is recommended to use uf_epoch for epoching, because you need to remove rows
% from EEG.unfold.X if the epoching function removed trials. Also cleaning
% of data is taken care of in uf_epoch
%
%Arguments:
%   cfg.method (string): (default pinv) 'glmnet','pinv','matlab','lsmr' are available. See
%      the uf_glmfit function for further information. By making use of
%      pinv, the linear model needs to be solved only once and can be
%      applied to all electrodes. The other solves iteratively solve for
%      each electtrode.
%
%  cfg.channel: (all) subselect a set of channels (in numbers, not strings)
%  cfg.glmnetalpha (1): used for glmnet --> see uf_glmfit
%  cfg.debug: used for matlab solver 
%  cfg.ica (boolean):0, use data or ICA components (have to be in
%               EEG.icaact). cfg.channel chooses the components.
%
%Return:
%  Returns a matrix (channel x pnts x predictors) of betas saved into
%  EEG.devon.beta
%
%*Example:*
% EEG = uf_glmfit_nodeconv(EEG)
%


cfg = finputcheck(varargin,...
    {'method',   'string', {'glmnet','pinv','matlab','lsmr'}, 'pinv';...
    'lsmriterations','integer',[],200;...
    'glmnetalpha','real',[],1;... # used for glmnet
    'debug','boolean',[],0;...
    'ica','boolean',[],0;...

    'channel','integer',[],1:size(EEG.data,1);
    },'mode','ignore');
if(ischar(cfg)); error(cfg);end

X = EEG.unfold.X;
if cfg.ica
    data = EEG.icaact;
else
    data = EEG.data;
end


%reshape instead of squeeze to keep the mxn matrix even when m or n is 1
if ~isfield(EEG.unfold,'timebasis')
    fprintf('deconvolution time-basis not found, assuming stick-functions / full \n')
    EEG.unfold.timebasis = eye(size(EEG.data,2));
end
    
beta = nan(size(data,1),size(data,2),size(X,2));
if strcmp(cfg.method,'pinv')
    %% Pseudoinverse
    % pseudoinverse solution can make use of basisFunctions. We therefore
    % overwrite the beta vector
    beta = nan(size(data,1),size(EEG.unfold.timebasis,1),size(X,2));

    Xinv = pinv(X);
    for c = cfg.channel
        % This is X^-1 * (time-basisFunction * Data) => We move the data first
        % in the same basis-domain i.e. splines, fourier or stick (which would
        % be the identity matrix) before running the regression
        beta(c,:,:)= (Xinv*(EEG.unfold.timebasis*reshape(data(c,:,:),[EEG.pnts, EEG.trials 1]))')';        
    end

    
    
elseif strcmp(cfg.method,'matlab') % save time
    warning('time-basis function currently not implemented')
    %% Matlab internal solver
    
    if cfg.debug
        spparms('spumoni',2)
    end
    for c = cfg.channel
        beta(:,:,c) = X \ squeeze(data(c,:,:))';
    end
    beta =     permute(beta,[3 2 1]);
    
elseif strcmp(cfg.method,'glmnet')
    warning('time-basis function currently not implemented')
    %% GLMNET
    beta = nan(size(EEG.data,1),size(data,2),size(X,2)+1); %plus one, because glmnet adds a intercept
    
    
    for e = cfg.channel
        t = tic;
        fprintf('\nsolving electrode %d (of %d electrodes in total)',e,length(cfg.channel))
        for time = 1:size(data,2)
            %glmnet needs double precision
            fit = cvglmnet(X,(double(squeeze(data(e,time,:))')),'gaussian',struct('alpha',cfg.glmnetalpha));
            
            %find best cv-lambda coefficients
            beta(e,time,:) = cvglmnetCoef(fit,'lambda_1se')';
            
            
            
        end
        fprintf('... took %.1fs',toc(t))
    end
    beta = beta([2:end 1],:,:); %put the dc-intercept last
    
    
elseif strcmp(cfg.method,'lsmr')
    warning('time-basis function currently not implemented')
    % go tru channels
    beta = nan(size(data,1),size(data,2),size(X,2)); %plus one, because glmnet adds a intercept
    
    for time = 1:size(data,2)
        for e = cfg.channel
            t = tic;
            fprintf('\nsolving electrode %d (of %d electrodes in total)',e,length(cfg.channel))
            
            % use iterative solver for least-squares problems (lsmr)
            data_tmp = squeeze(double(data(e,time,:)));
            [beta(e,time,:),ISTOP,ITN] = lsmr(X,data_tmp,[],10^-8,10^-8,[],cfg.lsmriterations); % ISTOP = reason why algorithm has terminated, ITN = iterations
            if ISTOP == 7
                warning(['The iterative least squares did not converge for channel ',num2str(e), ' after ' num2str(ITN) ' iterations'])
            end
            fprintf('... %i iterations, took %.1fs',ITN,toc(t))
            
            
        end
    end
    
    %% LSMR
    
end

EEG.unfold.beta_nodc = beta;
EEG.unfold.times = EEG.times/1000; %because seconds is better than ms!



end
