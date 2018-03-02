function [EEG] = uf_glmfit_nodc(EEG,varargin)
%A function to solve the inverse problem without deconvolution
% Currently it solves the pseudo-inverse of EEG.unfold.X and multiplies it
% to all channels.
%
%Arguments:
%   No arguments right now ('').
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
    
    'channel','integer',[],1:size(EEG.data,1);
    },'mode','ignore');
if(ischar(cfg)); error(cfg);end

X = EEG.unfold.X;


%reshape instead of squeeze to keep the mxn matrix even when m or n is 1
if ~isfield(EEG.unfold,'timebasis')
    warning('deconvolution time-basis not found, assuming stick-functions / full')
    EEG.unfold.timebasis = eye(size(EEG.data,2));
end
    
if strcmp(cfg.method,'pinv')
    %% Pseudoinverse
    Xinv = pinv(X);
    beta = calc_beta(EEG,Xinv);
    
elseif strcmp(cfg.method,'matlab') % save time
    warning('time-basis function currently not implemented')
    %% Matlab internal solver
    beta = nan(size(X,2),size(EEG.data,2),size(EEG.data,1)); %plus one, because glmnet adds a intercept
    if cfg.debug
        spparms('spumoni',2)
    end
    for c = cfg.channel
        beta(:,:,c) = X \ squeeze(EEG.data(c,:,:))';
    end
    beta =     permute(beta,[3 2 1]);
    
elseif strcmp(cfg.method,'glmnet')
    warning('time-basis function currently not implemented')
    %% GLMNET
    beta = nan(size(EEG.data,1),size(EEG.data,2),size(X,2)+1); %plus one, because glmnet adds a intercept
    
    
    for e = cfg.channel
        t = tic;
        fprintf('\nsolving electrode %d (of %d electrodes in total)',e,length(cfg.channel))
        for time = 1:size(EEG.data,2)
            %glmnet needs double precision
            fit = cvglmnet(X,(double(squeeze(EEG.data(e,time,:))')),'gaussian',struct('alpha',cfg.glmnetalpha));
            
            %find best cv-lambda coefficients
            beta(e,time,:) = cvglmnetCoef(fit,'lambda_1se')';
            
            
            
        end
        fprintf('... took %.1fs',toc(t))
    end
    beta = beta([2:end 1],:,:); %put the dc-intercept last
    
    
elseif strcmp(cfg.method,'lsmr')
    warning('time-basis function currently not implemented')
    % go tru channels
    beta = nan(size(EEG.data,1),size(EEG.data,2),size(X,2)); %plus one, because glmnet adds a intercept
    
    for time = 1:size(EEG.data,2)
        for e = cfg.channel
            t = tic;
            fprintf('\nsolving electrode %d (of %d electrodes in total)',e,length(cfg.channel))
            
            % use iterative solver for least-squares problems (lsmr)
            data = squeeze(double(EEG.data(e,time,:)));
            [beta(e,time,:),ISTOP,ITN] = lsmr(X,data,[],10^-8,10^-8,[],cfg.lsmriterations); % ISTOP = reason why algorithm has terminated, ITN = iterations
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

function [beta] = calc_beta(EEG,Xinv)

for c = 1:EEG.nbchan
    
    % This is X^-1 * (time-basisFunction * Data) => We move the data first
    % in the same basis-domain i.e. splines, fourier or stick (which would
    % be the identity matrix) before running the regression
    beta(c,:,:)= (Xinv*(EEG.unfold.timebasis*reshape(EEG.data(c,:,:),[EEG.pnts, EEG.trials 1]))')';

   
 

end
end