function r2 = uf_checkmodelfit(EEG,varargin)
% Function to compute R2, and partialR2 (in the future possibly AIC/BIC, but
% normally assumption is very likely wrong. Probably better to use CrossVal directly)
%
% Important: RESIDUALS ARE ONLY CALCULATED WHERE THERE ARE ENTRIES IN
% EEG.UNFOLD.XDC. thus _not_ in pauses/breaks/in-between trials etc.
%
% Input: EEG should the EEG either before or after uf_glmfit.
%
% Arguments:
% cfg.method (string): {Default: R2}
%                      -- R2: calculates coefficient of determination
%                      -- partialR2: calculates partial R2
%                          (model_with-model_without predictor) for all
%                           predictors (or if specified in cfg.variablename)
%                      -- crossValR2: R2 in cross-validated, folded by
%                           cfg.fold_event
%                      -- crossValpartialR2: same as partialR2 with
%                           CrossValidation
%
% cfg.fold_event (cell/string): The field that we can cut the data into
%                       folds without risking any overlap. E.g. this can be
%                       breaks, or if there is no overlap, trials.
%
% cfg.variablename (cell of strings): Specify the variablename
%                         (taken from EEG.unfold.variablenames) to which the
%                         partialR2 method should be applied to.
%
% cfg.channel:         On which channel(s) should the function be run?
%                      {Default: all channels}
%
% Related subfunctions: uf_cv_getfolds,

%% Input checks
cfg = finputcheck(varargin,...
    {'method', 'string',{'R2','partialR2','crossValR2','crossValpartialR2'}, 'R2';
    'fold_event','','',{}; %can be string or cell
    'channel','integer','',[];
    },'mode','ignore');
if(ischar(cfg)); error(cfg);end

% check EEG structure (continuous? required unfold fields?)
assert(ndims(EEG.data)==2,'EEG.data has to be continuous rather than epoched (i.e. 2 dimensions)')
assert(isfield(EEG,'unfold'),'Field EEG.unfold is missing')
assert(isfield(EEG.unfold,'Xdc'), 'Field EEG.unfold.Xdc is missing. Did you generate a (time-expanded) design matrix?')

% channel not specified: take all chans
if isempty(cfg.channel)
    cfg.channel = 1:size(EEG.data,1);
end

% debugging
% cfg.method     = 'crossValR2'
% cfg.method     = 'r2'
% cfg.fold_event = 'pause';

%% compute model fit (depending on method)
switch cfg.method
    
    case {'crossValR2','crossValpartialR2'}
        % get folds of the data for cross-validation
        % folds are defined based on the cfg.fold_event
        fprintf('\nPreparing data for crossvalidation...')
        [train,test] = uf_cv_getFolds(EEG,'fold_event',cfg.fold_event);
        r2 = []; % in case of crossValpartialR2 this is necessary
        
        fprintf('\nIdentified %i folds of the data (separated by event %s)',length(train),cfg.fold_event)
        % for each fold
        for fold = 1:length(train)
            fprintf('\nModeling fold: %i',fold)
            % select test & train time-ix
            
            % refit model for this fold
            EEGfold = EEG;
            EEGfold.unfold.Xdc = train(fold).Xdc;
            EEGfold = uf_glmfit(EEGfold);
            
            % R2 or partial R2?
            switch cfg.method
                case 'crossValpartialR2'
                    % potentially one could implement here a selection of
                    % only some variables. %% enhancement %%
                    r2_fold = partialR2(EEGfold,EEG.data(:,test(fold).ix),test(fold).Xdc(test(fold).ix,:),varargin);
                    r2_fold.fold = repmat(fold,size(r2_fold,1),1);
                    r2 = [r2; r2_fold];

                case 'crossValR2'
                    % calculate R2
                    r2(:,fold) = calc_r2(EEG.data(cfg.channel,test(fold).ix),test(fold).Xdc(test(fold).ix,:),EEGfold.unfold.beta_dc(cfg.channel,:,:));                   
            end
        end
        
    case 'partialR2'
        r2 = partialR2(EEG,EEG.data,EEG.unfold.Xdc,varargin);
        
    case 'R2'
        assert(isfield(EEG.unfold,'beta_dc'),'Missing field EEG.unfold.beta_dc. To obtain non-crossvalidated R2, please run uf_glmfit first manually')
        r2 = calc_r2(EEG.data(cfg.channel,:),EEG.unfold.Xdc,EEG.unfold.beta_dc(cfg.channel,:,:));
    
    otherwise
        error('Method not implemented, please check')
end

% debugging
% Specify / find Events that can be used for folding
% -- select one as a test-part
% -- temporarly blank the test-part of Xdc
% -- refit the model
% -- calculate R2
end

%% partialR2
% calculate partialR2 for variables (predictors)
function partial_r2_table = partialR2(EEG,testData,testXdc,input)

cfg = finputcheck(input,...
        {'variablename','cell','',{}; % cell of variablename names {'varnameA','varnameB'}
        'channel','integer',[],1:size(EEG.data,1);
        },'mode','ignore');
    if ischar(cfg); error(cfg);end

    % compute overall R2 of model
    r2_total = calc_r2(testData(cfg.channel,:),testXdc, EEG.unfold.beta_dc(cfg.channel,:,:)); % r2partial = r2_total - r2_without

    Xdc_terms2variablenames = EEG.unfold.cols2variablenames(EEG.unfold.Xdc_terms2cols); % which "Xdc" columns belong to which variables?
    varNames = []; partial_r2 = [];

    %% loop through each variable (predictor) in model
    for k = sort(unique(Xdc_terms2variablenames),'ascend')
        % in case we pre-specified variables to include, check whether variable is in the list, else skip
        if ~isempty(cfg.variablename)
            current_var = EEG.unfold.variablenames{k};

            if ~any(strcmp(current_var,cfg.variablename))
                fprintf('\nSkipping partialR2 computation for variablename: %s',current_var)
                continue
            end
            fprintf('\nCalculating partialR2 for variablename: %s',current_var)
        end
 
        % Now banish/delete the columns belonging to the variable from the 
        % design matrix and fit the model again
        % %% Enhancement %%: One could try permuting them similar to Mussal [et al] & Churchland 2019
        EEG_ca = EEG;
        
        ix_k = Xdc_terms2variablenames == k;
        EEG_ca.unfold.Xdc(:,ix_k) = [];
        EEG_ca.unfold.Xdc_terms2cols(ix_k) = [];
        EEG_ca.unfold.variablenames(k) = []; 
        
        cols_k = EEG_ca.unfold.cols2variablenames == k;
        EEG_ca.unfold.cols2eventtypes(cols_k) = [];
        EEG_ca.unfold.cols2variablenames(cols_k) = [];
        % fit model again
        EEG_ca = uf_glmfit(EEG_ca,'channel',cfg.channel); % debugging: betas of full model: 46*250*19

        % calculate partial R2
        testXdc_ca = testXdc;
        testXdc_ca(:,ix_k) = [];
        
        r2_ca = calc_r2(testData(cfg.channel,:),testXdc_ca,EEG_ca.unfold.beta_dc(cfg.channel,:,:)); % bug-fix by Olaf, 2020-08-29
        %r2_ca = calc_r2(testData(cfg.channel,:,:),testXdc_ca,EEG_ca.unfold.beta_dc(cfg.channel,:,:));
        partial_r2(:,end+1) = r2_total-r2_ca;
        varNames{end+1} = EEG.unfold.variablenames{k};
    end
    
    % collect partial R2 values
    partial_r2_table = [];
    for ch = 1:size(partial_r2,1)
        partial_r2_t          = table(repmat(cfg.channel(ch),size(partial_r2(ch,:)))',varNames',partial_r2(ch,:)','VariableNames',{'channel','variablenames','r2_ca'});
        partial_r2_t.r2_total = repmat(r2_total(ch),size(partial_r2_t,1),1);
        partial_r2_table = [partial_r2_table; partial_r2_t];
    end
end

%% calc_r2
% Calculate coefficient of determination for model
function r2 = calc_r2(data,Xdc,beta)
    model_ix  = any(Xdc,2); % rows with entries in DM (= all modeled samples/time points)
    residuals = (data' - Xdc*beta(:,:)')'; 
    residSS   = sum(residuals(:,model_ix).^2,2);
    totSS     = (sum(model_ix)-1) * var(data(:,model_ix),[],2);
    r2        = 1-residSS./totSS;
end