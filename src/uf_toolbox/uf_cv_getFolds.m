function [train,test] = uf_cv_getFolds(EEG,varargin)

cfg = finputcheck(varargin,...
    {'fold_event','','',{}; %can be string or cell
    'method','string',{'event','equalEntries'},'event';
    'nfold','integer',[],10; % used for equalEntries
    'allowOverlap','boolean',[],0; % only allows CV when there is no overlap between events (thus the folds are "independent")
    },'mode','ignore');

if cfg.method=="event"
    assert(~isempty(cfg.fold_event),'You have to specify a fold-event')
end

% find fold-breaks
switch cfg.method
    case "event"
        lat = {EEG.event.latency};
        if ischar(cfg.fold_event)
            cfg.fold_event ={cfg.fold_event};
        end
        lat = [lat{ismember({EEG.event.type},cfg.fold_event)}];
        
        % add first & last sample
        lat = [1 lat size(EEG.data,2)];
        lat = unique(lat); %remove double latencies
    case "equalEntries"
        ix = find(any(EEG.unfold.Xdc,2));
        lat = ix(round(linspace(1,length(ix),cfg.nfold+1)));
        
end

% check for no-overlap there
if ~cfg.allowOverlap
    assert(all(any(EEG.unfold.Xdc(round(lat),:),2)==0),'Error: The folding events have modelled overlap')
end


% for each fold
for fold = 1:(length(lat)-1)
    
    % select test & train time-ix
    test_ix = round(lat(fold)):round(lat(fold+1)-1);
    train_ix = [1:round(lat(fold))-1 round(lat(fold+1)):size(EEG.data,2)];
    % blank Xdc
    
    Xdc_test = EEG.unfold.Xdc;
    Xdc_test(train_ix,:) = 0;
    Xdc_train = EEG.unfold.Xdc;
    Xdc_train(test_ix,:) = 0;
    
    % check that there are any left in the fold
    
    try
        assert(any(any(Xdc_train,2)));
        assert(any(any(Xdc_test,2)));
    catch
        warning('found empty fold between latency %i and %i',round(lat(fold)),round(lat(fold+1)-1))
        continue
    end
    test(fold).Xdc = Xdc_test;
    train(fold).Xdc = Xdc_train;
    
    test(fold).ix = test_ix;
    train(fold).ix = train_ix;
    
end