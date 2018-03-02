function [EEG] = uf_imputeMissing(EEG,varargin)
% Deal with predictors for which some values are missing in design matrix
% You can either impute missing values or remove the predictors events for 
% which some values are missing
%
% Arguments:
% cfg.method:
%  * 'drop'     : (similar to R) Drop the whole event from the designmat (fill
%                 it with 0). This will lead to the event not being used
%                 for overlap correction!
%  * 'marginal' : fill in a random value from the marginal predictor-distribution
%                 in the future it might be interesting to implement not the
%                 marginal, but multivariate methods to conservate
%                 correlations between predictors (c.f. Horton & Kleinmann 2007)
%  * 'mean'     : fill in the mean value
%  * 'median'   : (Default) fill in the median value
%
% Returns:
% EEG.unfold.X in which missing NAN-values were imputed ('marginal', 'mean',
% 'median') or in which the events with missing predictor information were 
% removed ('drop'), which means put to 0

fprintf('\n%s(): Looking for missing predictor information...',mfilename)

cfg = finputcheck(varargin,...
    {'method',   'string', {'drop','marginal','mean','median'},'median';...
    },'mode','ignore');
if(ischar(cfg)); error(cfg);end

assert(length(EEG.event) == size(EEG.unfold.X,1),'error: EEG.event structure not as long as EEG.unfold.X')

%% find columns of design matrix with missing predictor information (NaNs)
missingColumn = find(any(isnan(EEG.unfold.X)));
if ~isempty(missingColumn)
    fprintf('\nMissing values in the following column(s): %s \n',strjoin(EEG.unfold.colnames(missingColumn),','))
else
    fprintf('\nNo missing values found in event structure\n')
    return
end

%% go tru columns containing missings
for pred = missingColumn
      
    %% find rows containing the event type
    cols2eventtypes = EEG.unfold.cols2eventtypes(pred);
    eventrows = strcmp(EEG.unfold.eventtypes{cols2eventtypes},{EEG.event(:).type})';
   
    %% get rows of this predictor with NaNs
    nanvals = isnan(EEG.unfold.X(:,pred));
    % consider only NaNs that are in "eventrows":
    nanIDX    =  nanvals & eventrows; 
    notnanIDX = ~nanvals & eventrows;
    
    %% get the "good" data (for interpolation)
    % take all values from predictor column that are in eventrows and *not* missing (not NaNs)
    X_pred = EEG.unfold.X(notnanIDX,pred);    
    
    nNAN = sum(nanIDX); % number of missing values
    perc = nNAN/length(eventrows)*100;  
    % bugfix Olaf: sum(eventrows) did not seem to make sense 
    % and gives the user a way too low feedback (e.g. 0.0001%) about the percentage of missing values in this predictor
    % eplaced by length(eventrows)

    %% user feedback
    if perc > 5
        warning([num2str(perc) '% of the values in predictor: ',EEG.unfold.colnames{pred},' are missing! This could bias your analysis'])
    end
    fprintf('\nimputing %.2f%% of values for predictor \"%s\" using method: \"%s\" \n',perc,EEG.unfold.colnames{pred},cfg.method) 
    % bugfix OD: this should not be an "else", since the feedback should also been given in the case above with > 5% missing
     
    %% now deal with missing values   
    switch cfg.method
        case 'drop'
            % in case of drop, we want to remove the whole event
            EEG.unfold.X(nanIDX,:) = 0; 
            continue
        case 'marginal'
            imputedValues(:) = datasample(X_pred,nNAN);
        case 'mean'
            imputedValues(:) = mean(X_pred);
        case 'median'
            imputedValues(:) = median(X_pred);
    end
    EEG.unfold.X(nanIDX,pred) = imputedValues(:);
end