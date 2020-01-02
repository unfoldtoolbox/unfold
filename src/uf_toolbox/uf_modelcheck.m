function r2 = uf_modelcheck(EEG,varargin)

cfg = finputcheck(varargin,...
    {'method', 'string',{'R2','commonalityR2','crossValR2','crossValcommonalityR2'}, 'R2';
    'fold_event','','',{}; %can be string or cell
    },'mode','ignore');
if(ischar(cfg)); error(cfg);end


assert(ndims(EEG.data)==2,'EEG.data has to be continuous, i.e. 2 dimensions')
assert(isfield(EEG,'unfold'))
assert(isfield(EEG.unfold,'Xdc'))

% cfg.method = 'crossValR2'
% cfg.method = 'r2'
% cfg.fold_event = 'pause';


switch cfg.method
    case {"crossValR2","crossValcommonalityR2"}
        
        [train,test]=uf_cv_getFolds(EEG,'fold_event',cfg.fold_event);
        r2 = []; % in case of crossValcommonalityR2 this is necessary
        
        % for each fold
        for fold = 1:length(train)
            fprintf('Running fold %i',fold)
            % select test & train time-ix
            
            % refit model
            EEGfold = EEG;
            EEGfold.unfold.Xdc = train(fold).Xdc;
            EEGfold = uf_glmfit(EEGfold);
            
            
            switch cfg.method
                case "crossValcommonalityR2"
                    % potentially one could implement here a selection of
                    % only some variables. %% enhancement %%
                    
                    r2_fold=  commonalityR2(EEGfold,EEG.data(:,test(fold).ix),test(fold).Xdc(test(fold).ix,:));
                    r2_fold.fold = repmat(fold,size(r2_fold,1),1);
                    r2 = [r2;r2_fold];
                    
                    
                case "crossValR2"
                    r2(fold) = calc_r2(EEG.data(:,test(fold).ix),test(fold).Xdc(test(fold).ix,:),EEGfold.unfold.beta_dc);
                    
                    %  calculate R2
                    
            end
            
        end
        
    case "commonalityR2"
        r2 = commonalityR2(EEG);
    case "R2"
        assert(isfield(EEG.unfold,'beta_dc'),'Please run uf_glmfit first manually to get non-crossvalidated R2')
        
        r2 = calc_r2(EEG.data,EEG.unfold.Xdc,EEG.unfold.beta_dc);
    otherwise
        error('not implemented method')
end



% Specify / find Events that can be used for folding
% select one as a test-part
% temporarly blank the test-part of Xdc
% refit the model
% calculate R2
end
function r2_commonality = commonalityR2(EEG,testData,testXdc)

if nargin == 3
    %crossvalidated r2
    r2_total    = calc_r2(testData,testXdc,      EEG.unfold.beta_dc);
    
else
    r2_total    = calc_r2(EEG.data,   EEG.unfold.Xdc,   EEG.unfold.beta_dc);
end
Xdc_terms2variablenames = EEG.unfold.cols2variablenames(EEG.unfold.Xdc_terms2cols); %which Xdc columns belong to which variables
varNames = [];
for k = sort(unique(Xdc_terms2variablenames),'ascend')
    EEG_ca = EEG;
    % banish these columns from the design matrix
    % Enhancement: One could try permuting them similar to Mussal [et al]...
    % & Churchland 2019. s
    ix_k = Xdc_terms2variablenames == k;
    EEG_ca.unfold.Xdc(:,ix_k) = [];
    EEG_ca.unfold.Xdc_terms2cols(ix_k) = [];
    %     EEG_ca.unfold.eventtypes(k) = [];
    EEG_ca.unfold.variablenames(k) = [];
    EEG_ca.unfold.cols2eventtypes(k) = [];
    EEG_ca.unfold.cols2variablenames(k) = [];
    
    EEG_ca = uf_glmfit(EEG_ca);
    
    if nargin == 3
        % crossvalidated R2
        testXdc_ca = testXdc;
        testXdc_ca(:,ix_k) = [];
        r2_ca = calc_r2(testData,testXdc_ca,EEG_ca.unfold.beta_dc);
        
        r2_commonality(k) = r2_total-r2_ca;
    else
        r2_ca = calc_r2(EEG.data,EEG_ca.unfold.Xdc,EEG_ca.unfold.beta_dc);
        
        r2_commonality(k) = r2_total-r2_ca ;
    end
    
    
    varNames{end+1} = EEG.unfold.variablenames{k};
    
end
r2_commonality = table(varNames',r2_commonality','VariableNames',{'variablenames','r2_ca'});
r2_commonality.r2_total = repmat(r2_total,size(r2_commonality,1),1);
end

% What has been modeled at all?
function r2 = calc_r2(data,Xdc,beta)
model_ix = any(Xdc,2);
residuals = (data' - Xdc*beta(:,:)')';

residSS = sum(residuals(:,model_ix).^2,2);
totSS = (sum(model_ix)-1) * var(data(:,model_ix),[],2);
r2 = 1-residSS./totSS;
end