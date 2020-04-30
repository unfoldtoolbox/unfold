function se = uf_se(EEG,varargin)
cfg = finputcheck(varargin,...
    {...
    'channels','real',1:size(EEG.data,1),[]; % combine these channels
    'restrictResidualToModelled','boolean',[],true; % calculate residuals only where something is modelled, or over all samples (including breaks etc.)
    'contrast','real','',[]; % whether to use raw data
    
    },'mode','error');
if ischar(cfg); error(cfg);end

% if we create the contrast, we can reshape the output to be
% size(beta_dc,[2,3])
reshapeToBetadc = false;
if isempty(cfg.contrast)
    % in case no custom contrast, we spit out the SE for each effect each
    % timepoint
    cfg.contrast = eye(size(EEG.unfold.Xdc,2));
    cfg.contrast = reshape(cfg.contrast,[],size(EEG.unfold.beta_dc,2),size(EEG.unfold.beta_dc,3));
    reshapeToBetadc = true;
end

% We can either calculate the residual variance based on ALL EEG sample
% (probably not what you want) or only on those part that were actually
% modelled (excluding breaks etc.)
if cfg.restrictResidualToModelled
    ix = any(EEG.unfold.Xdc,2);
else
    ix = 1:size(EEG.unfold.Xdc,1);
end
Xdc = EEG.unfold.Xdc(ix,:);

% Calculate residual variance
residualVar = var(mean(EEG.data(cfg.channels,ix),1)' - Xdc*mean(EEG.unfold.beta_dc(cfg.channels,:),1)');
warning('Autocorrelation was NOT taken into account. Therefore SE are UNRELIABLE. Use at your own discretion')

% Hat matrix
hat = full(inv(Xdc'*Xdc)) .* residualVar;

% apply contrast and convert from Variance to SD
se = sqrt(diag(cfg.contrast(:,:)*hat*cfg.contrast(:,:)'));

if reshapeToBetadc
    se = reshape(se,size(cfg.contrast,2),[]);
else
    % for custom contrast no guarantees can be made and we only reshape to
    % the number of contrasts - so we don't reshape
    
end