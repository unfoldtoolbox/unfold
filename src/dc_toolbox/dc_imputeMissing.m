function [EEG] = dc_imputeMissing(EEG,varargin)
% We impute missing values by several methods.
%
%Arguments:
% cfg.method:
%  * 'drop'     : (Default, similar to R) Drop the whole event from the designmat (fill
%               it with 0).
%  * 'marginal' : fill in a random value from the marginal predictor-distribution
%              in the future it might be interesting to implement not the
%              marginal, but multivariate methods to conservate
%              correlations between predictors (c.f. Horton & Kleinmann 2007)
%  * 'mean'     : fill in the mean value
%  * 'median'   : fill in the median value
%
%Returns:
% EEG.deconv.X is modified and the missing NAN-values were imputed

fprintf('\ndc_imputeMissing(): Handling missing predictor information...')

cfg = finputcheck(varargin,...
    {'method',   'string', {'drop','marginal','mean','median'},'drop';...
    },'mode','ignore');

if(ischar(cfg)); error(cfg);end


missingColumn = find(any(isnan(EEG.deconv.X)));

if ~isempty(missingColumn)
    fprintf('Missing values in the following column(s): %s \n',strjoin(EEG.deconv.colnames(missingColumn),','))
else
    display('No missing values found in event-structure')
end

for pred = missingColumn
    nanIDX = isnan(EEG.deconv.X(:,pred));
    cols2eventtypes = EEG.deconv.cols2eventtypes(pred);
    eventrows = strcmp(EEG.deconv.eventtypes{cols2eventtypes},{EEG.event(:).type}); % bugfix proposal by OD: field "EEG.deconv.event" does not exist
    X_pred = EEG.deconv.X(eventrows'&~nanIDX,pred);
    
    
    nNAN = sum(nanIDX);
    perc = nNAN/sum(eventrows)*100;
    if perc>5
        warning(['more than 5% of predictor (' num2str(perc) '%):',EEG.deconv.colnames{pred},' are missing. Thus could bias your analysis'])
    else
        fprintf('imputing %.2f%% of values for %s using method: %s \n',perc,EEG.deconv.colnames{pred},cfg.method)
    end
    
    imputedValues = nan(1,sum(nanIDX));
    switch cfg.method
        case 'drop'
            EEG.deconv.X(nanIDX,:) = 0;
            continue
        case 'marginal'
            imputedValues(:) = datasample(X_pred,nNAN);
        case 'mean'
            imputedValues(:) = mean(X_pred);
        case 'median'
            imputedValues(:) = median(X_pred);
    end
    EEG.deconv.X(nanIDX,pred) = imputedValues(:);
end