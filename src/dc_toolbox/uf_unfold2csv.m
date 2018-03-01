function varargout = dc_unfold2csv(unfold,varargin)
% Exports betas in an organized csv file to be opened in another tool
% returns a data-table
%
%Arguments:
%   cfg.deconv (boolean): Use the deconv betas (deconv.beta_dc) or the
%                         no-deconv betas(deconv.beta_nodc)
%   cfg.channel (integer): (Default: All channels) Limit to a list of specific channels
%
%   cfg.filename: filename for the csv file. if empty, only returns table
%
%Return:
%   Data-Table in the "tidy"-format:
%    Each observation (voltage/beta) has one row, channels, predictors etc.
%    gets one column


cfg = finputcheck(varargin,...
    { 'deconv', 'boolean',[],1;
    'channel','integer',[],[];
    'filename','string',[],'';
    },'mode','ignore');

if(ischar(cfg)); error(cfg);end

%% generate empty output structure

if cfg.deconv == 1
    data = unfold.beta;
elseif cfg.deconv == 0
    data = unfold.beta_nodc;
end

% remove channels with only nan (non-fitted)
rmchan = all(isnan(data(:,:)),2);

if ~isempty(cfg.channel)
    rmchan  = rmchan | ~ismember(1:size(data,1),cfg.channel);
end

fprintf('Exporting %i channel(s)\n',sum(~rmchan));
data = data(~rmchan,:,:);

nchan = size(data,1);
ntime = size(data,2);
npred = size(data,3);

predName= repmat({unfold.param(:).name}',1,ntime,nchan);
predName= permute(predName,[3 2 1]);

predValue= repmat([unfold.param(:).value]',1,ntime,nchan);
predValue= permute(predValue,[3 2 1]);

predEvent= repmat({unfold.param(:).event}',1,ntime,nchan);
predEvent= permute(predEvent,[3 2 1]);
predEvent = cellfun(@(x)strjoin(x),predEvent,'UniformOutput',0);


time   = repmat([unfold.times]',1,nchan,npred);
time   = permute(time,[2 1 3]);

if isfield(unfold,'chanlocs') && ~isempty(unfold.chanlocs)
    channels = repmat({unfold.chanlocs(~rmchan).labels}',1,ntime,npred);
else
    channels = repmat(1:nchan,1,ntime,npred);
end
% remove nans (due to unfitted data)


% construct table
t = table(predEvent(:),predName(:),predValue(:),channels(:),time(:),data(:),'VariableNames',{'event','predictor','predictorvalue','channel','time','data'});
if cfg.deconv == 1
    t.method = repmat({'deconvolution'},size(t,1),1);
elseif cfg.deconv == 0
    t.test = repmat({'no-deconvolution'},size(t,1),1);
end

%% Export
if ~isempty(cfg.filename)
    [path,name,ext] = fileparts(cfg.filename);
    if isempty(ext)
        ext = '.csv';
    end
    fullfilepath = fullfile(path,[name ext]);
    writetable(t,fullfilepath,'QuoteStrings',true)


    fprintf('sucessfully exported to %s\n',fullfilepath)
end

if nargout > 0
    varargout{1} = t;
end
