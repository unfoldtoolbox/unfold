
function [varargout] = uf_erpimage(EEG,varargin)

assert(isfield(EEG,'data'))
assert(ndims(EEG.data) == 2) % data have to be continuous
if 1 == 0
    %%
    cfgtmp = [];
    cfgtmp.type = 'deconv';
    cfgtmp.alignto = {'stimulus2'};
    %cfgtmp.remove = { {{'stimonset'},{}} };
    %cfgtmp.keep= { {{'stimulus2'},{'conditionA','2_(Intercept)'}} }; % possible: {'A','B'}, {{'A','B'}},{ {{'A'},{'B'}} , {{'C'},'D'}} }
    cfgtmp.remove= { {{'stimulus3'},{'2_(Intercept)'}} }; % possible: {'A','B'}, {{'A','B'}},{ {{'A'},{'B'}} , {{'C'},'D'}} }
    cfgtmp.sort_by = 'conditionA';
    cfgtmp.sort_align = 'stimulus2';
    cfgtmp.channel = 1;
    %cfgtmp.timelimits = [-0.1 0.5];
    cfgtmp.searchtime = [-10 10];
    
end
cfg = finputcheck(varargin,...
    {'type','string',{'raw','deconv','nodeconv','residuals','fullmodel','deconvDirect'},'deconv'; % deconv erpimage or not
    'remove','cell',[],{};   % which event(s) predictor(s) pair(s) to remove {'eventA',{'(Intercept)','stimA'}} (or a cell array of such cell arrays)
    'keep','cell',[],{};     % which event(s) predictor(s) pair(s) to keep   {'eventA',{'stimB'}} (or a cell array of such cell arrays)
    'alignto','cell',[],{};  % what should be the event to align the erpimage?
    'sort_align','',[],[];   % default same as alignto
    'sort_by','',[],'latency';   % could be in plotting function? What Eventfield to sort by
    'split_by','',[],[]; % make multiple subplots with separate ERPimages
    'winrej','real',[],[];
    'plot','boolean',[],0;  % if called without requesting output,
    'searchtime','real',[],[];
    'timelimits','real',[],[];
    'channel','integer',[],[];
    'figure','boolean',[],0;
    'caxis','real',[],[];
    },'mode','ignore');


if ischar(cfg); error(cfg);end
assert(isempty(cfg.keep)|isempty(cfg.remove),'either keep or remove, but dont specify both')
if isempty(cfg.sort_align)
    cfg.sort_align = cfg.alignto;
end
if cfg.figure
    figure
end
%%
if any(strcmp(cfg.type,{'raw','nodeconv','deconvDirect','residuals'}))
    fprintf('Epoching EEG in case of nodeconv')
    EEG_epoch = uf_epoch(EEG,struct('winrej',cfg.winrej,'timelimits',EEG.unfold.times([1,end])+[0 diff(EEG.unfold.times([1:2]))]));
    EEG_epoch = uf_glmfit_nodc(EEG_epoch,'channel',cfg.channel);
    ufresult= uf_condense(EEG_epoch);
else
    ufresult= uf_condense(EEG);
end

%%
events = [ufresult.unfold.eventtypes];
variable= [ufresult.unfold.variablenames]; % predictors


% default is to keep all
keep = 1:length(ufresult.unfold.colnames);


if ~(isempty(cfg.remove) && isempty(cfg.keep))
    % we fill it only with the keep variables
    list = [];
    
    if ~isempty(cfg.keep)
        fieldcontent = cfg.keep;
    else
        fieldcontent = cfg.remove;
    end
    for instance = fieldcontent
        current_event = instance{1}{1};
        ix_evt= find(cellfun(@(x)isequal(current_event,x),events));
        assert(length(ix_evt) == 1,'length of ix event was not 1')
        for field= instance{1}{2}
            ix_variable = find(cellfun(@(x)isequal(field{1},x),variable));
            assert(length(ix_variable) == 1,'length of ix variable was not 1')
            list = [list find(ufresult.unfold.cols2eventtypes == ix_evt & ufresult.unfold.cols2variablenames == ix_variable)];
        end
    end
    
    
    if ~isempty(cfg.keep)
        % keep the list only ("define" method)
        disp('Keeping the following predictors:')
        fprintf('%s, ',ufresult.unfold.colnames{list})
        fprintf('\n')
        keep = list;
    else
        % remove the list from all ("subtract" method)
        disp('Removing the following predictors:')
        fprintf('%s, ',ufresult.unfold.colnames{list})
        keep(list) = [];
        fprintf('\n')
    end
    
end

keep = sort(keep);
%%
switch cfg.type
    case 'raw'
        % nothing to do here
        %events, should be easy.
    case 'deconv'
        beta = ufresult.beta(cfg.channel,:,keep);
        assert(~all(isnan(beta(:))),'found only nans, did you specify & calculate the correct channel?')
        data = ufresult.unfold.Xdc(:,ismember(ufresult.unfold.Xdc_terms2cols,keep))*beta(:);
    case 'deconvDirect'
        beta = squeeze(ufresult.beta(cfg.channel,:,keep));
        data = ufresult.unfold.X(:,keep)*beta';
        data = data';
        
    case 'nodeconv'
        beta = squeeze(ufresult.beta_nodc(cfg.channel,:,keep));
        data = ufresult.unfold.X(:,keep)*beta';
        data = data';
    case {'residuals','fullmodel'}
        beta = ufresult.beta(cfg.channel,:,:);
        data = ufresult.unfold.Xdc*beta(:);
        
        
end

%%
% overwrite old data
EEG_new = EEG;
switch cfg.type
    case 'raw'
        
        EEG_modelled = EEG_epoch;
        EEG_modelled.data = EEG_modelled.data(cfg.channel,:,:);
    case {'nodeconv','deconvDirect'}
        EEG_modelled = EEG_epoch;
        EEG_modelled.data = permute(data,[3 1 2]);
    case {'deconv','residuals','fullmodel'}
        % first generate continuous data, for this remove the ICA info else
        % we get errors
        EEG_new.icaact = [];EEG_new.icasphere = [];EEG_new.icaweights= [];
        EEG_new.data = data';
        % now cut the data once more. This reintroduces overlapping potentials
        % - exactly what we want here
        EEG_modelled = uf_epoch(EEG_new,'timelimits',ufresult.unfold.times([1,end])+[0 diff(ufresult.times([1:2]))],'winrej',cfg.winrej);
end

if strcmp(cfg.type,'residuals')
    EEG_modelled.data = EEG_modelled.data -  EEG_epoch.data(cfg.channel,:,:);
end
[keep_epoch]= ~isnan(eeg_getepochevent(EEG_modelled,cfg.alignto,[0,0],'type')); % output in ms
EEG_modelled.data = EEG_modelled.data(:,:,keep_epoch);

sort_vector = eeg_getepochevent(EEG_modelled,cfg.sort_align,cfg.searchtime,cfg.sort_by); % output in ms
sort_vector = sort_vector(keep_epoch);

%% draw ERPimage
if isempty(cfg.caxis)
    cfg.caxis = prctile(EEG_modelled.data(:),[5 95]);
    cfg.caxis = [-max(abs(cfg.caxis)) max(abs(cfg.caxis))];
end
if nargout == 0 || cfg.plot == 1
    if isempty(cfg.split_by)
        
        
        erpimage(EEG_modelled.data,sort_vector,EEG_modelled.times,'',10,0,'caxis',cfg.caxis);
    else
        evt_tmp = {EEG.event.(cfg.split_by)};
        evt = evt_tmp(cellfun(@(x)~isempty(x),evt_tmp));
        splitlevel = unique(evt);
        n_splits = length(splitlevel);
        for n = 1:n_splits
            subplot(n_splits,1,n)
            evt_ix = strcmp(evt,splitlevel{n});
            erpimage(EEG_modelled.data(:,:,evt_ix),sort_vector(evt_ix),EEG_modelled.times,'',10,0,'caxis',cfg.caxis);
            title(sprintf('split by:%s, level:%s, effect of:%s',cfg.split_by,splitlevel{n},cfg.keep{1}{2}{:}))
        end
    end
else
    
    outdata = EEG_modelled.data;
    outsort = sort_vector;
    varargout{1} = outdata;
    varargout{2} = outsort;
end