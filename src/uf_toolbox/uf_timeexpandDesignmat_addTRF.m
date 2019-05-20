function EEG = uf_timeexpandDesignmat_addTRF(EEG,varargin)

%uf_designmat_addTRF
% Adds a timeexpanded continuous variable. The continuous variablce needs
% to be added as a new channel to the EEG.data matrix.
%
%Arguments:
%   timelimits (list of real): Similar to uf_timexpandDesignmat, how much
%              to timeexpand. This will be the resulting epoch size
%   channel (integer): Which continuous variable (aka channel) to use. If
%              channel is e.g. 65, the function expects the continuous variable to be
%              in EEG.data(65,:).
%   sparse (boolean): Should the resulting Xdc designmatrix be sparse
%   (recommended if you use other non-TRF predictors).
%   name (string): Name of the new predictor
%   
%
%Return:
%   EEG-Struct
%   * unfold.Xdc added EEG.data(cfg.channel,:) and timelags to the
%        designmatrix
%   * unfold.colnames (and others) added the cfg.name accordingly
% keep the designmatrix sparse? Put to 0

warning('EXPERIMENTAL FUNCTION - RISKY USAGE')
cfg = finputcheck(varargin,...
    {
    'timelimits','integer',[],[];...
    'channel', 'integer', [], [];...
    'sparse','boolean',[],1;... 
    'name','string',[],[];
    },'mode','ignore');

%XXX Check that the name is unique

if isfield(EEG,'unfold') && isfield(EEG.unfold,'times')
assert(min(EEG.unfold.times) == cfg.timelimits(1) &&max (EEG.unfold.times) == cfg.timelimits(2),'Timelimits do not match unfold.times min/max. currently only same timelimits allowed!')
end
assert(cfg.channel <= size(EEG.data,1),'Chosen channel is larger than size(EEG.data,1)')

if isfield(EEG,'unfold') && isfield(EEG.unfold,'X')
    error('Please run uf_timeexpandDesignmat before running uf_timeexpandDesignmat_addTRF')
end
    
if ~isfield(EEG,'unfold') || ~isfield(EEG.unfold,'Xdc')
    fprintf(['Could not find Xdc, assuming only TRFs are fitted and initializing the Xdc field \n'...
    'Important: If you want to combine rERPs and TRFs you HAVE to use uf_timeexpand BEFORE uf_timeexpandDesignmat_addTRF\n'])

%     fnlist = {'formula','eventtypes','cols2eventtypes','variablenames','variabletypes','splines','colnames','X','cols2variablenames'};
    EEG.unfold.Xdc = nan(size(EEG.data,2),0);
    EEG.unfold.colnames = {};
    EEG.unfold.Xdc_terms2cols = [];
    EEG.unfold.eventtypes = {};
    EEG.unfold.cols2eventtypes = [];
end
tmin = cfg.timelimits(1)*EEG.srate;
tmax = cfg.timelimits(2)*EEG.srate;
additionalColumns = lagGen(EEG.data(cfg.channel,:)',tmin:tmax); % lagGentaken from mTRF toolbox
% XXX Combine lagGen with 
% additionalColumns = lagGen(rand(size(EEG.data,2),1),tmin:tmax); % taken from mTRF toolbox
for c = additionalColumns
    EEG = uf_designmat_addcol(EEG,c,cfg.name);
end
if ~isfield(EEG.unfold,'times')
    EEG.unfold.times =  (tmin:tmax)/EEG.srate;
end