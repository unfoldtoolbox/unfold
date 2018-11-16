function [EEG, eventwise] = uf_continuousArtifactExclude(EEG,varargin)
%% Function to exclude (artifactual) continuous data from being modeled
% This function expects a rejection vector and excludes those intervals from
% being modeled in the design matrix. That means all predictor values in 
% the given intervals are set to 0 in the time-expanded design matrix and 
% therefore ignored in the modeling process
%
%Arguments:
%   cfg.winrej (integer): A (2xn) array with n from-to pairs of samples to
%   be excluded from further processing This is the same output as from
%   EEGlabs' eegplot rej
%
%Return:
%   EEG-Structure
%   * unfold.X: All elements between the from-to pairs got set to 0

%   eventwise
%   structured array containing infos about how much data was excluded

%Example:
% We want to exclude three intervals (10:50, 100:120 etc.) of data that 
% contain artifacts (marked manually or by some automatic function like 
% uf_continuousArtifactDetect() ):
%| cfgReject = [];
%| cfgReject.winrej = [10,50; 100,120; 300,310];
%| EEG = uf_artefactRemoveDesignmat(cfgReject,EEG)


cfg = finputcheck(varargin,...
    {'winrej',   'integer', [], [];...
    },'mode','error');
if(ischar(cfg)); error(cfg);end

rej = [];

cfg.winrej = round(cfg.winrej);

times = EEG.unfold.times;

eventlat = [EEG.event.latency];
for k = 1:size(cfg.winrej,1)
    rej = [rej cfg.winrej(k,1):cfg.winrej(k,2)];
end

% structure that saves how many trials were excluded
eventwise = struct('type',unique({EEG.event.type}));

fprintf('Portions of data removed split up by each event in EEG.event\n')
% add how much total time the events occupy (in s)
for ev = 1:length(eventwise)
    evIDX = strcmp({EEG.event.type},eventwise(ev).type);
    
    eventwise(ev).unit = 'seconds';
    eventwise(ev).total_signallength = size(EEG.data,2) * 1/EEG.srate;
    
    % do as if no overlap
    eventwise(ev).total_evttime= sum(evIDX) * range(EEG.unfold.times);
    
    % include the overlap
    activesamp = zeros(size(EEG.data,2),1);
    for evlat = [EEG.event(evIDX).latency]
        ix = round((evlat+min(EEG.unfold.times)*EEG.srate ) : (evlat+max(EEG.unfold.times)*EEG.srate));
        ix(ix<=0) = []; %negative time does not exist
        ix(ix>size(EEG.data,2)) = []; %larger than EEG does not exist
        activesamp( ix ) = 1;
    end
    eventwise(ev).abs_evttime= sum(activesamp) * 1/EEG.srate;
    
    rejsamp = zeros(size(EEG.data,2),1);
    for k = 1:size(cfg.winrej,1)
        rejwindow =  cfg.winrej(k,1):cfg.winrej(k,2);
        rejsamp(rejwindow) = 1;
    end
    eventwise(ev).abs_removed = sum(rejsamp&activesamp) * 1/EEG.srate;
    
    fprintf('Type: %-14s modelled eventtime: %9.2fs \t rejected eventtime: %9.2fs \t percent removed: %4.1f%%\n',eventwise(ev).type,eventwise(ev).abs_evttime,eventwise(ev).abs_removed, eventwise(ev).abs_removed / eventwise(ev).abs_evttime * 100)
end

EEG.unfold.Xdc(round(rej),:) = 0;
fprintf('\nRemoving %.2f%% of rows from design matrix (filling them with zeros) \n',length(unique(rej))/size(EEG.unfold.Xdc,1)*100)

end
