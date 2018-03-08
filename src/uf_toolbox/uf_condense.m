function output = uf_condense(EEG,varargin)
%% Condense results in new structure. Apply timebasis (if necessary) 
%Returns an "ufresult"-struct that contains the predictor betas over time
%and accompanying information. This structure is further used in all
%plotting functions. This function also applys the timebasis (if you
%specified something else than the default 'stick' in
%uf_timeexpandDesignmat()
%
%Arguments:
%   EEG(struct): A Struct containing EEG.unfold.beta_dc
%   cfg.deconv (integer): 1, use EEG.unfold.beta_dc, the deconvolved betas
%                         0, use EEG.unfold.beta_nodc, betas without
%                         deconvolution
%                         -1 (default), autocheck which fields are avaiable
%                         and returns both
%   cfg.channel(array): Restrict the beta-output to a subset of
%                         channels. Default is all channels
%
%Return:
%   ufresult.beta= (nchans x time x parameters)
%   ufresult.beta_nodc = (nchans x time x parameters) (only if unfold=0 or -1)
%   ufresult.param = (struct size: parameters) each field contains the values of the respective parameter.
%   ufresult.unfold = EEG.unfold
%   ufresult.times = EEG.times
%   ufresult.chanlocs = EEG.chanlocs
%
%**Example:**
%
%ufresult = uf_condense(EEG)
%
%ufresult.param(X):
%
%* name: name of the variable, e.g.: 'continuousA'
%* value: value of the predictor, e.g. '50'
%* event: event of the variable, e.g.: 'eventA'

assert(isfield(EEG.unfold,'beta_nodc')|isfield(EEG.unfold,'beta_dc'),'Input Error: Could not find beta-estimates. Did you run uf_glmfit?')
if isfield(EEG.unfold,'beta_nodc')
    nchan = size(EEG.unfold.beta_nodc,1);
elseif isfield(EEG.unfold,'beta_dc')
    nchan = size(EEG.unfold.beta_dc,1);
end
    
cfg = finputcheck(varargin,...
    { 'deconv', 'integer',[-1 0 1],-1;
      'channel','integer',[],1:nchan;
    },'mode','error');

if(ischar(cfg)); error(cfg);end

if isfield(EEG,'nbchan')
    assert(all(ismember(cfg.channel,1:nchan)))
end

beta_dcExists   = isfield(EEG.unfold,'beta_dc')&& isnumeric(EEG.unfold.beta_dc);
beta_nodcExists = isfield(EEG.unfold,'beta_nodc') && isnumeric(EEG.unfold.beta_nodc);

if cfg.deconv == 1
    assert(beta_dcExists,'beta_dc missing or not numeric')
elseif cfg.deconv== 0
    assert(beta_nodcExists,'beta_nodc missing or not numeric')
elseif cfg.deconv == -1 % auto detect, recursive call
    assert(beta_dcExists | beta_nodcExists,'either beta_dc or beta_nodc need to exist. Did you fit the model already?')
    %-------------- Recursive part
    if beta_dcExists && beta_nodcExists
        output1 = uf_condense(EEG,'channel',cfg.channel,'deconv',1);
        output2 = uf_condense(EEG,'channel',cfg.channel,'deconv',0);

        output = output1;
        output.beta_nodc = output2.beta_nodc;
    elseif beta_dcExists
        output = uf_condense(EEG,'channel',cfg.channel,'deconv',1);
    else
        output = uf_condense(EEG,'channel',cfg.channel,'deconv',0);
    end
    return
    %-------------- End Recursive part
end
%% generate empty output structure


% find out which (if any) parameters are modeled by splines, we only want to add the
% combined spline estimate, not each individually.

[splIdxList,paramList] = uf_getSplineidx(EEG);


paramList = 1:size(EEG.unfold.X,2);



signal = nan(length(cfg.channel),length(EEG.unfold.times),length(paramList));
times = EEG.unfold.times;

% initialize vectors
value = nan(1,size(signal,3));
name = cell(1,size(signal,3));
event = name;
type = name;

loopRunner = 1;
%% go trough all parameters in parameterList
for pred = paramList
    if cfg.deconv
        signal(:,:,loopRunner) = EEG.unfold.beta_dc(cfg.channel,:,pred)*EEG.unfold.timebasis;
    else
        signal(:,:,loopRunner) = EEG.unfold.beta_nodc(cfg.channel,:,pred)*pinv(EEG.unfold.timebasis)';
    end
    
    % change name incase of spline and no conversion
    if ismember(pred,splIdxList)
        colname = EEG.unfold.colnames(pred);
        
        
        splt = regexp(colname,'([-]?[\d]*\.?[\d]*)$','tokens');
        number = splt{1}{1}{1};
        name(loopRunner) = {colname{1}(1:(end-length(number)-1))}; % -1 to get rigd of the "_"
        value(loopRunner) = str2num(number);
        %             name(loopRunner)= splt(1);
    else
        value(loopRunner) = nan;
        name(loopRunner)= EEG.unfold.colnames(pred);
    end
    event(loopRunner)= EEG.unfold.eventtypes(EEG.unfold.cols2eventtypes(pred));
    type(loopRunner) = EEG.unfold.variabletypes(EEG.unfold.cols2variablenames(pred));
    loopRunner = loopRunner+1;
%     end

end

output = struct();
output.unfold = EEG.unfold;
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
