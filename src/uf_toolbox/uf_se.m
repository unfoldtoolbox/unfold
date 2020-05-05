function [se,contrast] = uf_se(EEG,varargin)
cfg = finputcheck(varargin,...
    {...
    'channels','real',[1,size(EEG.data,1)],[]; % combine these channels
    'restrictResidualToModelled','boolean',[],true; % calculate residuals only where something is modelled, or over all samples (including breaks etc.)
    'contrast','real','',[]; % whether to use raw data
    'spline_addmarginal','boolean','',0; % whether to add marginal of spline effects
    
    },'mode','error');
if ischar(cfg); error(cfg);end

if isempty(cfg.channels)
    warning('No channels selected. Calculating SE over all channels.Are you sure this is what you want?')
    cfg.channels = 1:size(EEG.data,1);
end
% if we create the contrast, we can reshape the output to be
% size(beta_dc,[2,3])
reshapeToBetadc = false;
if isempty(cfg.contrast)
    % in case no custom contrast, we spit out the SE for each effect each
    % timepoint
    cfg.contrast = eye(size(EEG.unfold.Xdc,2));
    cfg.contrast = reshape(cfg.contrast,[],size(EEG.unfold.beta_dc,2),size(EEG.unfold.beta_dc,3));
    reshapeToBetadc = true;
    assert(~cfg.spline_addmarginal,'you have to specify a contrast matrix if you want to add spline marginals')
else
    assert(size(cfg.contrast,2) == size(EEG.unfold.beta_dc,2))
    assert(size(cfg.contrast,3) == size(EEG.unfold.beta_dc,3))
    
    if cfg.spline_addmarginal
        % Often when we have splines, we want to add those back into our
        % SE-calculation.
        % This is only supported if not all contrasts want to be
        % calculated. Then it would be a bit messy, I'm sure
        warning('Functionallity alpha-version. Adding add-marginal to all non-zero contrast values (per contrast)!')
        % loop over all specified contrasts
        for c = 1:size(cfg.contrast,1)
            % for each contrast loop over all splines
            for spl = EEG.unfold.splines
                % what values do we have to set the splines, so that we add
                % the mean(splineValue) to all contrasts?
                tmp = spl{1}.splinefunction(mean(spl{1}.paramValues),spl{1}.knots);
                tmp(spl{1}.removedSplineIdx) = []; % remove the removed spline
                
                % figure out to which columns this spline needs to be added
                var_ix = find(strcmp(spl{1}.name,EEG.unfold.variablenames));
                splIx= EEG.unfold.cols2variablenames == var_ix;
                
                % figure out at which timepoints we have to add the
                % marginal. This is where cfg.contrast is not 0.
                event_ix = unique(EEG.unfold.cols2eventtypes(splIx));
                assert(length(event_ix)==1,'Error, found multiple events where I should have found only one. Please report this at https://github.com/unfoldtoolbox/unfold/issues')
                time_ix = any(cfg.contrast(c,:,EEG.unfold.cols2eventtypes==event_ix),3);
                
                % add the mean-spline designmatrix entry to the right place
                cfg.contrast(c,time_ix,splIx) = tmp;
            end
        end
        
    end
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
assert(~isnan(residualVar),'residual Variance was NaN')
% Hat matrix
tic
XtX = Xdc'*Xdc;
toc
hat = full(inv(XtX)) .* residualVar;
toc

% apply contrast and convert from Variance to SD
se = sqrt(diag(cfg.contrast(:,:)*hat*cfg.contrast(:,:)'));

if reshapeToBetadc
    se = reshape(se,size(cfg.contrast,2),[]);
else
    % for custom contrast no guarantees can be made and we only reshape to
    % the number of contrasts - so we don't reshape
    
end
contrast = cfg.contrast;  % to get the currently used contrast back out