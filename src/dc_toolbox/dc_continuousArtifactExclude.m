function EEG = dc_continuousArtifactExclude(EEG,varargin)
%Function to exclude (artifactual) continuous data from being modeled
% This function inputs a rejection vector and excludes the content from
% being modeled in the design matrix. That means it sets all predictor
% values at the given times to 0.
%
%Arguments:
%   cfg.winrej (integer): A (2xn) array with n from-to pairs of samples to be excluded from further processing
%
%Return:
%   EEG-Structure
%
%   * deconv.X: All elements between the from-to pairs got set to 0
%
%Example:
% We want to exclude three sections that are supposedly artifactual
%| cfgReject = [];
%| cfgReject.winrej = [10,50; 100,120; 300,310];
%| EEG = dc_artefactRemoveDesignmat(cfgReject,EEG)

cfg = finputcheck(varargin,...
    {'winrej',   'integer', [], [];...
    'zerodata','boolean', [],0;... % undocumented, also removes portions in EEG.data. This is useful sometimes because e.g. the stopping-criterion of the LSMR iterative solver depends on the data. With huge outliers, strange criterions can exist
    },'mode','ignore');
if(ischar(cfg)); error(cfg);end

rej = [];
cfg.winrej = round(cfg.winrej);
for k = 1:size(cfg.winrej,1)
    rej = [rej cfg.winrej(k,1):cfg.winrej(k,2)];
end

fprintf('\nremoving %.1f%% from design matrix (fill it with zeros) \n',length(unique(rej))/size(EEG.deconv.dcX,1)*100)
EEG.deconv.dcX(round(rej),:) = 0;

if cfg.zerodata
EEG.data(:,round(rej)) = 0;
end

end
