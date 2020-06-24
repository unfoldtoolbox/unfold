function test_continuousArtifact()
%% Check continuousArtifactDetect

EEG = simulate_test_case(5,'noise',2,'basis','box');
EEG = uf_designmat(EEG,'formula','y~1+continuousA','eventtypes','stimulusA');

cfgClean = [];
cfgClean.amplitudeThreshold=150;
cfgClean.windowsize=1000; %in ms
cfgClean.stepsize=100; %in ms
cfgClean.combineSegments=200; %in ms
cfgClean.channels=1;

EEG.data(2,:) = EEG.data(1,:);
EEG.data(1,500) = 1000;
EEG.data(1,505) = 1000;
EEG.data(1,525) = 1000;
EEG.data(2,525) = 1000;
EEG.data(2,625) = 1000;
EEG = eeg_checkset(EEG);
winrej = uf_continuousArtifactDetect(EEG,cfgClean);

resultMat = [492 513; 517 533];

if ~near(winrej,resultMat)
    error('problems with continuous artefact detection')
end
EEG = uf_timeexpandDesignmat(EEG,'timelimits',[-1,2]);

%% Check continuousArtifactExclude
% now we delete it based on the continous design matrix
EEG2 = uf_continuousArtifactExclude(EEG,'winrej',winrej);

t1 = ~all(all(EEG2.unfold.Xdc(winrej(1,1):winrej(1,2),:)==0));
t2 = ~all(all(EEG2.unfold.Xdc(winrej(2,1):winrej(2,2),:)==0));
if t1 || t2
    error('problems with continuous artefact rejection')
end

%% Test probabilistic exclusion of artefacts
EEG = simulate_test_case(5,'srate',100,'noise',2,'basis','box');

EEG.data(2,:) = EEG.data(1,:)+5*randn(size(EEG.data));

rng(1)
for k  = 1:5000:size(EEG.data,2)
    EEG.data(1:2,k:k+300) = 100+randn(2,301)*10;
end


EEG = eeg_checkset(EEG);

winrej_loc = uf_continuousJointProbArtifactDetect(EEG,'locthresh',3,'globthresh',3000);
winrej_glob = uf_continuousJointProbArtifactDetect(EEG,'locthresh',30000,'globthresh',3);

winrej= uf_continuousJointProbArtifactDetect(EEG,'locthresh',3,'globthresh',3);

% combine the two separate measures
winrej_combine = uf_combineWinrej(winrej_loc,winrej_glob);
% fuse together the winrej' that were calculated concurrently
winrej_fuse = uf_combineWinrej(winrej);
% check that they are equal
assert(all(winrej_combine(:) == winrej_fuse(:)))
%% Some tests to check if optional arguments work
winrej = uf_continuousJointProbArtifactDetect(EEG,'robust_normalization',0);
winrej = uf_continuousJointProbArtifactDetect(EEG,'channel',2);
winrej = uf_continuousJointProbArtifactDetect(EEG,'verbose',1);
%% Visual Check
% Visual check that the functions identify the massively obvious artefacts,
% additional parts can potentially be identified as well
figure,

for k = 1:size(winrej,1)
    
    h = fill(winrej(k,[1 1 2 2]),[-200 300 300 -200],'black');
    set(h,'facealpha',0.3)
    hold on
end
for k = 1:size(winrej_loc,1)
    h = fill(winrej_loc(k,[1 1 2 2]),[-200 300 300 -200],'red');
    set(h,'facealpha',0.3)
end
for k = 1:size(winrej_glob,1)
    h = fill(winrej_glob(k,[1 1 2 2]),[-200 300 300 -200],'blue');
    set(h,'facealpha',0.3)
end

plot(EEG.data')

%% Check uf_combineWinrej

% check if it removes overlap
A = [1 5;
    3 8;
    10 15;
    15 20;
    25 30];
X = uf_combineWinrej(A);
assert(all(X(:) == [1 10 25 8 20 30]'))

% check unsorted
A = [3 8;
    1 5;
    15 20;
    25 30
    10 15;];
X = uf_combineWinrej(A);
assert(all(X(:) == [1 10 25 8 20 30]'))

% Combine
A = [1 5;
    3 8;];
B = [10 15;];
C = [15 20;
    25 30];
X = uf_combineWinrej(A,B,C);
assert(all(X(:) == [1 10 25 8 20 30]'))