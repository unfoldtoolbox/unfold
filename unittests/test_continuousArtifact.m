function test_continuousArtifact()
%%
EEG = simulate_test_case(5,'noise',0,'basis','box');
EEG = uf_designmat(EEG,'formula','y~1+continuousA','eventtypes','stimulusA');

cfgClean = [];
cfgClean.amplitudeThreshold=150;
cfgClean.windowsize=1000; %in ms
cfgClean.stepsize=100; %in ms
cfgClean.combineSegements=200; %in ms
cfgClean.channels=1;


EEG.data(1,500) = 1000;
EEG.data(1,505) = 1000;
EEG.data(1,525) = 1000;
winrej = uf_continuousArtifactDetect(EEG,cfgClean);

resultMat = [492 513; 517 533];

if ~near(winrej,resultMat)
    error('problems with continuous artefact detection')
end
EEG = uf_timeexpandDesignmat(EEG,'timelimits',[-1,2]);

% now we delete it based on the continous design matrix
EEG2 = uf_continuousArtifactExclude(EEG,'winrej',winrej);

t1 = ~all(all(EEG2.deconv.Xdc(winrej(1,1):winrej(1,2),:)==0));
t2 = ~all(all(EEG2.deconv.Xdc(winrej(2,1):winrej(2,2),:)==0));
if t1 || t2
    error('problems with continuous artefact rejection')
end
