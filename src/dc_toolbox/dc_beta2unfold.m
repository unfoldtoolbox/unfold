function output = dc_beta2unfold(EEG,varargin)
%Returns an "unfold"-struct that contains the predictor betas over time
%and accompanying information. This structure is further used in all
%plotting functions
%
%Arguments:
%   EEG(struct): A Struct containing EEG.deconv.beta_dc
%   cfg.deconv (integer): 1, use EEG.deconv.beta_dc, the deconvolved betas
%                         0, use EEG.deconv.beta_nodc, betas without
%                         deconvolution
%                         -1 (default), autocheck which fields are avaiable
%                         and returns both
%   cfg.channel(array): Restrict the beta-output to a subset of
%                         channels. Default is all channels
%   cfg.pred_value(cell): Only necessary if splines are used. One entry per parameter:
%       {{'par1',[10 20 30]},{'par2',[0,1,2]}}.
%       This evaluates parameter 1 at the values 10,20 and 30. Parameter 2
%       at 0, 1 and 2. Default behaviour: evaluates 7 linearly spaced values between the min + max. of the
%       parameterdomain
%
%Return:
%   unfold.beta= (nchans x time x parameters)
%   unfold.beta_nodc = (nchans x time x parameters) (only if deconv=0 or -1)
%   unfold.param = (struct size: parameters) each field contains the values of the respective parameter.
%   unfold.deconv = EEG.deconv
%   unfold.times = EEG.times
%   unfold.chanlocs = EEG.chanlocs
%
%**Example:**
%
%unfold = dc_beta2unfold(EEG)
%
%unfold.param(X):
%
%* name: name of the variable, e.g.: 'continuousA'
%* value: value of the predictor, e.g. '50'
% * event: event of the variable, e.g.: 'eventA'


if isfield(EEG.deconv,'beta_nodc')
    nchan = size(EEG.deconv.beta_nodc,1);
elseif isfield(EEG.deconv,'beta_dc')
    nchan = size(EEG.deconv.beta_dc,1);
end
    
cfg = finputcheck(varargin,...
    { 'deconv', 'integer',[-1 0 1],-1;
      'channel','integer',[],1:nchan;
      'convertSplines','',[],[];
    },'mode','error');

if(ischar(cfg)); error(cfg);end

if ~isempty(cfg.convertSplines)
    error('convertSplines has been deprecated. See issue #17')
end
if isfield(EEG,'nbchan')
    assert(all(ismember(cfg.channel,1:nchan)))
end

beta_dcExists   = isfield(EEG.deconv,'beta_dc')&& isnumeric(EEG.deconv.beta_dc);
beta_nodcExists = isfield(EEG.deconv,'beta_nodc') && isnumeric(EEG.deconv.beta_nodc);

if cfg.deconv == 1
    assert(beta_dcExists,'beta_dc missing or not numeric')
elseif cfg.deconv == 0
    assert(beta_nodcExists,'beta_nodc missing or not numeric')
elseif cfg.deconv == -1 % auto detect, recursive call
    assert(beta_dcExists | beta_nodcExists,'either beta_dc or beta_nodc need to exist. Did you fit the model already?')
    %-------------- Recursive part
    if beta_dcExists && beta_nodcExists
        output1 = dc_beta2unfold(EEG,'channel',cfg.channel,'deconv',1);
        output2 = dc_beta2unfold(EEG,'channel',cfg.channel,'deconv',0);

        output = output1;
        output.beta_nodc = output2.beta_nodc;
    elseif beta_dcExists
        output = dc_beta2unfold(EEG,'channel',cfg.channel,'deconv',1);
    else
        output = dc_beta2unfold(EEG,'channel',cfg.channel,'deconv',0);
    end
    return
    %-------------- End Recursive part
end
%% generate empty output structure


% find out which (if any) parameters are modeled by splines, we only want to add the
% combined spline estimate, not each individually.

[splIdxList,paramList] = dc_getSplineidx(EEG);
% nSpl = sum(ismember(paramList,splIdxList));
% if nSpl == 0
%     nSplineBetas = 0;
% else
%     nSplineBetas = -nSpl + sum( cellfun(@(x)size(x.basis,1),EEG.deconv.predictorSplines) );
% end
% if ~cfg.convertSplines
paramList = 1:size(EEG.deconv.X,2);
% end


signal = nan(length(cfg.channel),length(EEG.deconv.times),length(paramList));
times = EEG.deconv.times;

% initialize vectors
value = nan(1,size(signal,3));
name = cell(1,size(signal,3));
event = name;
type = name;

loopRunner = 1;
%% go trough all parameters in parameterList
for pred = paramList
    if cfg.deconv
        signal(:,:,loopRunner) = EEG.deconv.beta_dc(cfg.channel,:,pred)*EEG.deconv.timebasis;
    else
        signal(:,:,loopRunner) = EEG.deconv.beta_nodc(cfg.channel,:,pred)*pinv(EEG.deconv.timebasis)';
    end
    
    % change name incase of spline and no conversion
    if ismember(pred,splIdxList)
        colname = EEG.deconv.colnames(pred);
        
        
        splt = regexp(colname,'([-]?[\d]*\.?[\d]*)$','tokens');
        number = splt{1}{1}{1};
        name(loopRunner) = {colname{1}(1:(end-length(number)-1))}; % -1 to get rigd of the "_"
        value(loopRunner) = str2num(number);
        %             name(loopRunner)= splt(1);
    else
        value(loopRunner) = nan;
        name(loopRunner)= EEG.deconv.colnames(pred);
    end
    event(loopRunner)= EEG.deconv.eventtype(EEG.deconv.cols2eventtype(pred));
    type(loopRunner) = EEG.deconv.variableType(EEG.deconv.cols2variableNames(pred));
    loopRunner = loopRunner+1;
%     end

end

output = struct();
output.deconv = EEG.deconv;
if cfg.deconv
    output.beta = signal;
else
    output.beta_nodc = signal;
end
output.times = times;
if isfield(EEG,'chanlocs')
    output.chanlocs = EEG.chanlocs;
else
    warning('no chanlocs found')
end
output.param = struct('value',num2cell(value),'name',name,'event',event,'type',type);
if isfield(EEG,'chanlocs')&&~isempty(EEG.chanlocs)
    output.chanlocs= EEG.chanlocs(cfg.channel);
end
