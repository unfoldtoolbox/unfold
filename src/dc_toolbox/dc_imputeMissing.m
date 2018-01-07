function [EEG] = dc_imputeMissing(EEG,varargin)
% Deal with predictors for which some values are missing in design matrix
% You can either impute missing values or remove the predictors events for 
% which some values are missing
%
% Arguments:
% cfg.method:
%  * 'drop'     : (Default, similar to R) Drop the whole event from the designmat (fill
%                 it with 0).
%  * 'marginal' : fill in a random value from the marginal predictor-distribution
%                 in the future it might be interesting to implement not the
%                 marginal, but multivariate methods to conservate
%                 correlations between predictors (c.f. Horton & Kleinmann 2007)
%  * 'mean'     : fill in the mean value
%  * 'median'   : fill in the median value
%
% Returns:
% EEG.deconv.X in which missing NAN-values were imputed ('marginal', 'mean',
% 'median') or in which the events with missing predictor information were 
% removed ('drop')

fprintf('\n%s(): Looking for missing predictor information...',mfilename)

cfg = finputcheck(varargin,...
    {'method',   'string', {'drop','marginal','mean','median'},'drop';...
    },'mode','ignore');
if(ischar(cfg)); error(cfg);end

assert(length(EEG.event) == size(EEG.deconv.X,1),'error: EEG.event structure not as long as EEG.deconv.X')

%% find columns of design matrix with missing predictor information (NaNs)
missingColumn = find(any(isnan(EEG.deconv.X)));
if ~isempty(missingColumn)
    fprintf('\nMissing values in the following column(s): %s \n',strjoin(EEG.deconv.colnames(missingColumn),','))
else
    fprintf('\nNo missing values found in event structure\n')
    return
end

%% go tru columns containing missings
for pred = missingColumn
      
    %% find rows containing the event type
    cols2eventtypes = EEG.deconv.cols2eventtypes(pred);
    eventrows = strcmp(EEG.deconv.eventtypes{cols2eventtypes},{EEG.event(:).type})';
   
    %% get rows of this predictor with NaNs
    nanvals = isnan(EEG.deconv.X(:,pred));
    % consider only NaNs that are in "eventrows":
    nanIDX    =  nanvals & eventrows; 
    notnanIDX = ~nanvals & eventrows;
    
    %% get the "good" data (for interpolation)
    % take all values from predictor column that are in eventrows and *not* missing (not NaNs)
    X_pred = EEG.deconv.X(notnanIDX,pred);    
    
    nNAN = sum(nanIDX); % number of missing values
    perc = nNAN/length(eventrows)*100;  
    % bugfix Olaf: sum(eventrows) did not seem to make sense 
    % and gives the user a way too low feedback (e.g. 0.0001%) about the percentage of missing values in this predictor
    % eplaced by length(eventrows)

    %% user feedback
    if perc > 5
        warning([num2str(perc) '% of the values in predictor: ',EEG.deconv.colnames{pred},' are missing! This could bias your analysis'])
    end
    fprintf('\nimputing %.2f%% of values for predictor \"%s\" using method: \"%s\" \n',perc,EEG.deconv.colnames{pred},cfg.method) 
    % bugfix OD: this should not be an "else", since the feedback should also been given in the case above with > 5% missing
     
    %% now deal with missing values   
    switch cfg.method
        case 'drop'
            EEG.deconv.X(nanIDX,:) = 0; 
            % bugfix olaf: this index was previously wrong, since nanIDX was created only on the "eventrows" subset of EEG.deconv.X, not the whole EEG.deconv.X!
            continue % shouldn't this be "return"?
        case 'marginal'
            imputedValues(:) = datasample(X_pred,nNAN);
        case 'mean'
            imputedValues(:) = mean(X_pred);
        case 'median'
            imputedValues(:) = median(X_pred);
    end
    EEG.deconv.X(nanIDX,pred) = imputedValues(:);
end