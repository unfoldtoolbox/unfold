function EEG_epoch = dc_epoch(EEG,varargin)
%% Epoch the data according to the deconv structure
% Deconvolution works on continuous data, thus to compare it to the
% "normal" use-case, we have to epoch it. Because the data has not been
% cleaned yet, we do this in this function. We additionally remove trials
% from deconv.X that were removed during epoching.
%
%Arguments:
%   cfg.winrej (integer): A (2xn) array with n from-to pairs of samples to be excluded from further processing
%   cfg.timelimits (float): min+max of the epoch in seconds
%   EEG (eeglab): the EEG set, need to have EEG.deconv.Xdc compatible with the size of EEG.data
%
%Returns:
%   Epoched EEG file to cfg.timelimits
%
%*Example:*
% EEG_epoch = dc_epoch(EEG,'winrej',winrej,'timelimits',cfgTimeexpand.timelimits)
%


cfg = finputcheck(varargin,...
    { 'winrej', 'integer',{}, [];...
    'timelimits', 'integer',{}, [];...
    },'mode','ignore');
if(ischar(cfg)), error(cfg); end

%% Remove events in windows with bad (continuous) EEG data
% This routine checks whether an event is contained in the rejection window
% If yes, it removes it from further analyses
if ~isempty(cfg.winrej)
    removeEvent = zeros(1,length(EEG.event));
   
    % go tru "bad data" windows
    for w = 1:size(cfg.winrej,1)
        win = cfg.winrej(w,:); % window (in samples)
        eL = round([EEG.event(:).latency]); % event latencies (in samples)
        % go tru events
        for e = 1:length(eL)
            if removeEvent(e)==1
                continue % we know already it is overlapping
            end
            eW = eL(e)+(cfg.timelimits*EEG.srate);

            if win(1) < eW(1)
                removeEvent(e) = win(2) >= eW(1);
            else
                removeEvent(e) = eW(2) >= win(1);
            end
        end
    end
    
    
    fprintf('Deleting %i events due to data cleaning \n',sum(removeEvent))
    % delete the event, but also the designmat entry
    EEG.event(removeEvent==1) = [];
    EEG.deconv.X(removeEvent==1,:) = [];
end

%% Epoch the data
EEG.urevent = EEG.event;
[EEG_epoch,event_ind] = pop_epoch(EEG,[EEG.deconv.eventtypes{:}],cfg.timelimits);
fprintf('Recalculating the EEG_epoch field using eeg_checkset, might take some time \n')
EEG_epoch             = eeg_checkset(EEG_epoch,'eventconsistency'); % needed to add EEG_epoch.epoch field, don't ask why



% we need to remove the entries in the table that got removed due to the
% epoching step (e.g. boundary events). But we only have ind, which is relative to the eventtypes.
% Thus we need to get the eventtypes once more and select only the matching
% ones

eventType = [EEG.deconv.eventtypes{:}];
eventType = eventType(~cellfun(@(x)isnan(x(1)),eventType));
convertIND = find(ismember({EEG.event(:).type},eventType));
EEG_epoch.deconv.X = EEG_epoch.deconv.X(convertIND(event_ind),:);

EEG_epoch.urevent = EEG.urevent(event_ind);