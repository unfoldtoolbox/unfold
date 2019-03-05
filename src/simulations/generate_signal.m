function EEG=generate_signal(X,eventTimes,sig,signals,cfg)

if any(eventTimes<= 0|eventTimes >= cfg.pnts)
    warning('some eventtimes are smaller than 0 or longer than the data length')
end

% Generate the signal

separateSignal = zeros(length(signals),length(sig.shape));
sigCount = 0;

data =  zeros(cfg.pnts,1);
for runner = 1:length(eventTimes)
    signalidx = eventTimes(runner):eventTimes(runner)+length(sig.shape)-1;
    cutsignalIX = signalidx<cfg.pnts;
    
    loopSig = sig.shape(cutsignalIX);
    loopIdx = signalidx(cutsignalIX);
    currSig = zeros(length(signals),length(loopSig));
    for s = 1:length(signals)
        
        switch signals(s).type
            case 'intercept'
                assert(X(runner,s)==1) % data.X(s,runner) has to be 1 for the intercept anyway
                currSig(s,:) = loopSig .* signals(s).effectsize .* X(runner,s); 
            case {'1x2','continuous'}
                currSig(s,:) = loopSig .* signals(s).effectsize .* X(runner,s); 
            case 'spline'
                currSig(s,:) = loopSig .* signals(s).effectsize .* signals(s).function(X(runner,s)); 
        end
    end
        if ~isempty(currSig)
            sigCount = sigCount + 1;
            separateSignal(:,cutsignalIX) = separateSignal(:,cutsignalIX) + (currSig-separateSignal(:,cutsignalIX))/sigCount;  % incremental mean
        end
        data(loopIdx)=  data(loopIdx)+sum(currSig,1)';
    
end

if cfg.noise
    data = data + cfg.noise*randn(size(data));
end

% make pseudo EEG-lab interface
EEG = eeg_emptyset;


EEG.sim.separateSignal = separateSignal;
EEG.sim.noOverlapSignal = sum(EEG.sim.separateSignal,1);

EEG.data = permute(data,[3 1 2]);% we want "electrode (1)" x "epochtime" x "epochs"
EEG.pnts = size(EEG.data,2);
EEG.trials = size(EEG.data,3);
EEG.nbchan = size(EEG.data,1);
EEG.srate = cfg.srate;

EEG.sim.signals = signals;
EEG.sim.sig =sig;
EEG.sim.X = X;
EEG.sim.eventTimes= eventTimes;

EEG.times = 1/EEG.srate:1/EEG.srate:EEG.pnts/EEG.srate;

EEG.event = struct('latency',mat2cell(eventTimes',ones(length(eventTimes),1)), ...
    'type',signals(1).eventname);


for e = 1:length(EEG.event)
    for t =1:size(X,2)
        EEG.event(e).(signals(t).predictorName) = X(e,t);
    end
    
end
