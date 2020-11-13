function [EEGcont,mask] = uf_epoch2cont(EEG)
if 1 == 0
    
   EEG = pop_epoch(a.EEG,[],[-0.1,1]); 
   EEGcont = uf_epoch2cont(EEG);
end
% function to concatenate epochs back into a continuous datastream
% requires correct EEG.urevent. Use at own precaution!
warning('No unit tests yet, please write them before usage')

%We need to know how much to shift the signal "backwards"
emin = sum(EEG.times<0);
elength = length(EEG.times) ;

% Find the latencies of the epoching event
epochUrevent = eeg_getepochevent(EEG,'',[0 0],'urevent');
epochingsample = [EEG.urevent(epochUrevent).latency];

%contdata should range from
tmin = 1;
% latency of last epoch + epoch length
tmax = elength+ epochingsample(end);
contdata = nan(size(EEG.data,1),tmax);
% loop over events and fill in data
for e = 1:length(EEG.epoch)
    estart = epochingsample(e) - emin;
    ix = estart:estart+elength-1;
    contdata(:,ix) = EEG.data(:,:,e);
end
mask = isnan(contdata(1,:));
EEGcont = eeg_emptyset();
EEGcont.data = contdata;
EEGcont.event = EEG.urevent;
EEGcont  = eeg_checkset(EEGcont);