function [varargout] = uf_erpimage(EEG,varargin)
% This function produces ERPimages with/without modelled data.
% DESCRIPTION
%
%Arguments:
% Mandatory
%  cfg.channel (integer):   Which channel(s) should the erpimage be plotted of?
%
% Specify to-be-plotted data
%  cfg.method(string): default 'deconv'.
%
%        * 'raw'        : Using raw-data
%        * 'modelled'   : Using modelled data (y_hat), usually deconv or
%                         no-deconv version (specified in cfg.datafield)
%        * 'residual'   : Plot residuals (useful to check for non-modelled
%                         stimulus locked activity)
%
%  cfg.datafield (string):   default 'beta'. Any field of "EEG.unfold",
%                            typically beta or beta_nodc - only in custom
%                            cases something else
%
%  cfg.alignto (String/cell of strings): Which event should mark t=0 in the
%                                        erpimage? (don't confuse with
%                                        sort_alignto). Default is all
%                                        events specified in unfold.X
%  cfg.overlap (boolean): default 0. Re-introduce the overlap? This is
%                         useful to check the modelfit, or useful in cases
%                         where cfg.keep/cfg.remove is used.
%
%  cfg.addResiduals (integer 0-2):  default 0. 
%                                   If 1: Adds the overlap-including residuals.
%                                         That is: y_cont - X_dc*beta
%                                   If 2: Adds the corresponding residuals,
%                                         this will be identical to "cfg.method='raw'"
%                                         except when cfg.keep/cfg.remove
%                                         is used
%
%  cfg.winrej (2 integer):   default []. rejection matrix used to fit the
%                           model. necessary if you want to remove the same noisy trials you removed
%                           during modelfitting
%
% Modify y_hat / modelled data
%
%  cfg.keep (cell):   default {}. Requires a: {{{},{}}} (yes, a cell of cell of cells - sorry :|
%                           Possibilty to specify what event-predictor combinations to keep.
%                           In case of "event1:~1+condA, event2:~1+condB"
%                           {{'event1',{'condA'}},{'event2',{'2_(Intercept)'}}}
%                           would plot single trials without
%                           event1:(Intercept) and event2:condB
%
%  cfg.remove (cell):   default {}. Similar to keep, but specify only the
%                           combinations you want to remove.
%
%
% Sorting the erpimage
%  cfg.sort_by ():   default ''. XXXmeter is not used.
%
%  cfg.sort_alignto (string/cell of strings):   default 'cfg.alignto'. Event(s) the
%                               ERPimage should be aligned to (specifies 
%                               where to look for information for the
%                               black-sorting line). By default alignto
%                               aligns to the erpimage-event (cfg.alignto)
%
%  cfg.sort_by (string):   default 'latency'. Looks in EEG.event.(cfg.by) for
%                           the value to sort by. Has to be an eventname of
%                           EEG.event
%
%  cfg.sort_time (2 integer):   default 'whole epoch'. You can subset where to
%                     look for "cfg.sort_alignto" events. Three typical
%                     use-cases: 
%                     1) cfg.time[0 1.5], in case the sort_alignto event 
%                        of interest occurs after the cfg.alignto event
%                     2) [-0.3 0] in case the event is before (check out cfg.sort_direction in this case)
%                     3) cfg.time=[0,0], in case the erpimage should be sorted by a field in cfg.alignto
%
%  cfg.sort_direction (string):   default 'forward'. Should the first or
%             last event in cfg.sort_time be selected, in case there are multiple
%             ones? Can be 'forward' or 'backward'. Useful if you are looking for the
%             first cfg.sort_alignto event prior to the cfg.alignto event
%
%
% Plotting options
%  cfg.plot (2 integer):  Force the plot. By default the plot is surpressed
%                       if you want the output
%  cfg.sort_time (2 integer):   default 'whole epoch'. You can subset where to
%  cfg.split_by (string):   default []. Use subplots to split the
%                       erpimage by a categorical variable. Directly uses
%                       EEG.event.(cfg.split_by)
%  cfg.figure (integer):   default: 0, plot into a new figure?
%  cfg.caxis (2 integer):   Specify caxis, by default let eeglab decide

%
%Returns:
%   * (optional) data: The erpimage data, ntimes x ntrials
%   * (optional) sort: the sorting index used to sort the erpimage
%
%*Example:*
%

assert(isfield(EEG,'data'),'uf_erpimage needs the EEG file after uf_glmfit, before uf_condense')
% assert(ismatrix(EEG.data)) % data have to be continuous

cfg = finputcheck(varargin,...
    {...%#'method','string',{'raw','deconv','nodeconv','residuals','fullmodel','deconvDirect'},'deconv'; %
    'type','string',{'raw','modelled','residual'},'modelled'; % whether to use raw data
    'datafield','string',[],'beta'; % which field to plot, popular: beta & beta_nodc
    'overlap','boolean',[],0;   % re-introduce the overlap
    'addResiduals','integer',[0:2],0; % add residuals?
    ...
    ... %in case of "modelled
    'remove','cell',[],{};     % which event(s) predictor(s) pair(s) to remove {'eventA',{'(Intercept)','stimA'}} (or a cell array of such cell arrays)
    'keep','cell',[],{};       % which event(s) predictor(s) pair(s) to keep   {'eventA',{'stimB'}} (or a cell array of such cell arrays)
    ...
    ... % what to cut & sort to
    'alignto','',[],{};    % what should be the event to align the erpimage?
    'sort_alignto','',[],[];     % default same as alignto (i.e. if you want to search for alignmentevents starting from a different event)
    'sort_by','',[],'latency'; % could be in plotting function? What Eventfield to sort by
    'sort_time','real',[],[];  % when to when to look for the sort_align event
    'sort_direction','string',{'forward','backward'},'forward'; % in case of multiple events, take first or last?
    ...
    ... % Plot Options
    'split_by','',[],[];       % make multiple subplots with separate ERPimages
    'winrej','real',[],[];     % remove parts of the continuous data
    'plot','boolean',[],nargout==0;     % if called without requesting output,
    'timelimits','real',[],[]; %from when to when
    'channel','integer',1:size(EEG.data,1),[];
    'figure','boolean',[],0;
    'caxis','real',[],[];
    },'mode','error');


if ischar(cfg); error(cfg);end
assert(isempty(cfg.keep)|isempty(cfg.remove),'either keep or remove, but dont specify both')
assert(~isempty(cfg.channel),'please specify a channel (or set of channels) for the ERPimage')

if isstr(cfg.alignto)
    cfg.alignto = {cfg.alignto};
end
assert(iscell(cfg.alignto))
if isempty(cfg.sort_alignto)
    cfg.sort_alignto = cfg.alignto;
end
if cfg.figure
    figure
end
%%
if ismatrix(EEG.data)
    % we should epoch the data
    [EEG_epoch] = uf_epoch(EEG,struct('winrej',cfg.winrej,'timelimits',EEG.unfold.times([1,end])+[0 diff(EEG.unfold.times([1:2]))]));
    EEG_epoch = uf_glmfit_nodc(EEG_epoch,'channel',cfg.channel);
    ufresult= uf_condense(EEG_epoch);
else
    EEG_epoch = EEG;
    ufresult= uf_condense(EEG_epoch);
    
end
assert(isfield(ufresult,cfg.datafield),sprintf('"datafield":%s, not found',cfg.datafield))
%% Find which parameters to keep in the model
events  = [ufresult.unfold.eventtypes];
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
        current_event = instance{1}(1);
        ix_evt= find(cellfun(@(x)isequal(current_event,x),events));
                    

        assert(length(ix_evt) == 1,'Could not find event:%s, did you misspell something?',current_event{1})
        for field= instance{1}{2}
            ix_variable = find(cellfun(@(x)isequal(field{1},x),variable));
            assert(length(ix_variable) == 1,'Could not find event:%s, field:%s, did you misspell something?',current_event{1},field{1})
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
if cfg.type == "modelled" || cfg.type == "residual"
    beta = ufresult.(cfg.datafield)(cfg.channel,:,keep);
    beta = squeeze(mean(beta,1));
    if cfg.type == "residual" || cfg.addResiduals
        beta_full = ufresult.(cfg.datafield)(cfg.channel,:,:);
        beta_full = squeeze(mean(beta_full,1));
    end
    assert(~all(isnan(beta(:))),'found only nans, did you specify & calculate the correct channel?')
    
    if cfg.overlap
        % go back to the continuous domain and multiply via Xdc
        data = ufresult.unfold.Xdc(:,ismember(ufresult.unfold.Xdc_terms2cols,keep))*beta(:); % time x 1
        if cfg.type == "residual" || cfg.addResiduals
            data_yhat = ufresult.unfold.Xdc*beta_full(:);
        end
    else
        if length(keep) == 1
            beta = beta'; % because of the "nice" behaviour of matlab to turn every Nx1 to a 1xN we have to undo that
        end
        % directly use unfold.X
        data = ufresult.unfold.X(:,keep)*beta';
        data = data'; % time x trial
         if cfg.type == "residual" || cfg.addResiduals
            
            data_yhat = ufresult.unfold.X*beta_full';
            data_yhat = data_yhat';
         end
    end
    
    
    % add a first (singleton) dimension mimicking a single channel
    data = permute(data,[3 1 2]);
    if exist("data_yhat","var")
        data_yhat = permute(data_yhat,[3 1 2]);
    end

    
elseif cfg.type == "raw"
    data = mean(EEG_epoch.data(cfg.channel,:,:),1);
end
%% Generate a EEG structure with the new data,
EEG_new = eeg_emptyset;
EEG_new.srate = EEG.srate;
EEG_new.unfold = EEG.unfold;
EEG_new.event = EEG.event;

if ~cfg.overlap || cfg.type == "raw"
    % in case of no overlap (or raw), the new data are already epoched
    EEG_out = EEG_epoch;
    EEG_out.data = data;
    keepIX = 1:size(EEG_out.data,3); 
else
    % first generate continuous data
    EEG_new.data = data;
    % now cut the data once more. This reintroduces overlapping potentials
    % - exactly what we want here
    [EEG_out,keepIX] = uf_epoch(EEG_new,'timelimits',ufresult.unfold.times([1,end])+[0 diff(ufresult.times([1:2]))],'winrej',cfg.winrej);
    
end
%%
% Adding back residuals
if cfg.addResiduals || cfg.type == "residual"
    assert(cfg.type~="raw",'cannot add residuals in case of "raw"')
    fprintf('adding residuals\n')
    EEG_residuals = EEG_new;
    % this is the fully modelled data
    EEG_residuals.data  = data_yhat;
    if cfg.overlap
        %in case of overlap we have continuous residuals and have to cut
        %them once more to epochs
        EEG_new.event = EEG_epoch.event;
        EEG_residuals = uf_epoch(EEG_residuals,'timelimits',ufresult.unfold.times([1,end])+[0 diff(ufresult.times([1:2]))],'winrej',cfg.winrej);
    end
    
    % calculate residuals as y-yhat    
    EEG_residuals.data = mean(EEG_epoch.data(cfg.channel,:,:),1) - EEG_residuals.data;

    %%
    if cfg.type == "residual"
        % show only residuals
        EEG_out.data = EEG_residuals.data;
    else
        % add residuals back in
        EEG_out.data = EEG_out.data + EEG_residuals.data;
    end
    clear EEG_residuals % free RAM
end


%% check out which epoch to keep
[keep_epoch]= ~isnan(eeg_getepochevent(EEG_out,cfg.alignto,[0,0],'type')); % output in ms
assert(sum(keep_epoch)>0,'Did not find any events to align to in the epochs')
EEG_out.data = EEG_out.data(:,:,keep_epoch);
fprintf('Aligning erpimage to event %s, %i epochs found\n',strjoin(cfg.alignto,':'),sum(keep_epoch))

[sort_vector,sort_vector_cell] = eeg_getepochevent(EEG_out,cfg.sort_alignto,cfg.sort_time*1000,cfg.sort_by); % output in ms


sort_isempty = cellfun(@(x)isempty(x),sort_vector_cell);
switch cfg.sort_direction
    case 'forward'
        sort_tmp = cellfun(@(x)x(1),sort_vector_cell(~sort_isempty));
    case 'backward'
        sort_tmp = cellfun(@(x)x(end),sort_vector_cell(~sort_isempty));
    otherwise
        error('unspecified sortdirection')
end
fprintf('Sorting erpimage from event %s %ss by field %s\n',strjoin(cfg.sort_alignto,':'),cfg.sort_direction,cfg.sort_by)

sort_vector = nan(size(sort_vector));% not strictly necessary, but maybe I made a mistake somewhere
sort_vector(~sort_isempty) = sort_tmp;
% select only those that we want to keep anyway

sort_vector = sort_vector(keep_epoch);

%% draw ERPimage
if isempty(cfg.caxis)
    cfg.caxis = prctile(EEG_out.data(:),[5 95]);
    cfg.caxis = [-max(abs(cfg.caxis)) max(abs(cfg.caxis))];
end

if cfg.plot == 1
    if isempty(cfg.split_by)
        
        
        erpimage(EEG_out.data,sort_vector,EEG_out.times,'',10,0,'caxis',cfg.caxis);
    else
        evt_tmp = {EEG.event.(cfg.split_by)};
        evt = evt_tmp(cellfun(@(x)~isempty(x),evt_tmp));
        splitlevel = unique(evt);
        n_splits = length(splitlevel);
        for n = 1:n_splits
            subplot(n_splits,1,n)
            evt_ix = strcmp(evt,splitlevel{n});
            erpimage(EEG_out.data(:,:,evt_ix),sort_vector(evt_ix),EEG_out.times,'',10,0,'caxis',cfg.caxis);
            title(sprintf('split by:%s, level:%s, effect of:%s',cfg.split_by,splitlevel{n},cfg.keep{1}{2}{:}))
        end
    end
    cmap = cbrewer('div','RdBu',256);
    colormap(cmap(end:-1:1,:));
elseif nargout>0
    
    outdata = EEG_out.data;
    outsort = sort_vector;
    varargout{1} = outdata;
    varargout{2} = outsort;
end
