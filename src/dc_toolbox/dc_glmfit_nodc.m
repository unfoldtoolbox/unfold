function [EEG] = dc_glmfit_nodc(EEG,varargin)
%A function to solve the inverse problem without deconvolution
% Currently it solves the pseudo-inverse of EEG.deconv.X and multiplies it
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
% EEG = dc_glmfit_nodeconv(EEG)
%


cfg = finputcheck(varargin,...
    {'method',   'string', {'glmnet','pinv','matlab'}, 'pinv';...
    'glmnetalpha','real',[],1;... # used for glmnet
    'debug','boolean',[],0;...

    'channel','integer',[],1:size(EEG.data,1);
    },'mode','ignore');
if(ischar(cfg)); error(cfg);end

X = EEG.deconv.X;

if strcmp(cfg.method,'pinv')
    
    Xinv = pinv(X);
    beta = calc_beta(EEG,Xinv);
    
elseif strcmp(cfg.method,'matlab') % save time
    
    beta = nan(size(X,2),size(EEG.data,2),size(EEG.data,1)); %plus one, because glmnet adds a intercept
    if cfg.debug
        spparms('spumoni',2)
    end
    for c = cfg.channel
     beta(:,:,c) = X \ squeeze(EEG.data(c,:,:))';
    end
    beta =     permute(beta,[3 2 1]);
    
elseif strcmp(cfg.method,'glmnet')
    
    beta = nan(size(X,2)+1,size(EEG.data,2),size(EEG.data,3)); %plus one, because glmnet adds a intercept
    for time = 1:size(EEG.data,2)
        for e = cfg.channel
            t = tic;
            fprintf('\nsolving electrode %d (of %d electrodes in total)',e,length(cfg.channel))
            %glmnet needs double precision
             fit = cvglmnet(X,(double(squeeze(EEG.data(e,time,:))')),'gaussian',struct('alpha',cfg.glmnetalpha));
            
            %find best cv-lambda coefficients
            beta(:,e,time) = cvglmnetCoef(fit,'lambda_1se')';
            
            fprintf('... took %.1fs',toc(t))
            
        end
    end
    beta = beta([2:end 1],:,:); %put the dc-intercept last
    
    
    
    
end

EEG.deconv.XBeta = beta;
EEG.deconv.dcBasistime = EEG.times/1000; %because seconds is better than ms!



end

function [beta] = calc_beta(EEG,Xinv)

% beta = nan(EEG.nbchan,size(EEG.deconv.dcBasis,1),size(Xinv,1));
% beta = nan(EEG.nbchan,size(EEG.data,2),size(Xinv,1));
% warning('dc_glmfit_nodc: deactivated basis functions for nodc-glmfit for now.')
for c = 1:EEG.nbchan
    
    % This is X^-1 * (time-basisFunction * Data) => We move the data first
    % in the same basis-domain i.e. splines, fourier or stick (which would
    % be the identity matrix) before running the regression
    
    %reshape instead of squeeze to keep the mxn matrix even when m or n is 1
    if ~isfield(EEG.deconv,'dcBasis')
        warning('deconvolution time-basis not found, assuming stick-functions / full')
        beta(c,:,:)= (Xinv*(reshape(EEG.data(c,:,:),[EEG.pnts, EEG.trials 1]))')';
    else
        
        beta(c,:,:)= (Xinv*(EEG.deconv.dcBasis*reshape(EEG.data(c,:,:),[EEG.pnts, EEG.trials 1]))')';
    end
    %     beta(c,:,:)= (Xinv*(reshape(EEG.data(c,:,:),[EEG.pnts, EEG.trials 1]))')';
end

end