function test_glmfit
%%
sprev = rng(1);

% Basis functions
intercept = struct();
intercept.eventname = 'stimulusA';
intercept.type = 'intercept';
intercept.predictorName = 'intercept';
intercept.overlap = 0; % note that 0, means that the AVERAGE overlap is exactly one signal-length (i.e. 1s in our simulations). That means if we do not want ANY overlap, we specify -1 here!
intercept.effectsize = 3;

cont = struct();
cont.eventname = 'stimulusA';
cont.type = 'continuous';
cont.overlap = 0;
cont.predictorName = 'continuousA';
cont.effectsize = 1;
cont.range = [0,100];

signals{1} = intercept;
signals{1}.range = nan;
signals{1}.overlap = -1;
signals{1}(2) = cont;
signals{1}(2).overlap = -1;
signals{1}(2).effectsize= 2.5;
signals{1}(3) = cont;
signals{1}(3).predictorName = 'continuousB';
signals{1}(3).overlap = -1;
signals{1}(3).effectsize= -1.5;

EEG = simulate_data(signals,'noise',0);

rng(sprev);

%%

cfgDesign = [];

cfgDesign.coding = 'dummy';
cfgDesign.formula   = {'y~1+continuousA+continuousB'};
cfgDesign.eventtype = {'stimulusA'};

EEG = dc_designmat(EEG,cfgDesign);
% add a second channel
% (this only works because there is no overlap + a single event)
EEG.data(2,:) = (EEG.data);
EEG = eeg_checkset(EEG);



EEG = dc_timeexpandDesignmat(EEG,'timelimits',[-1 2]);
%%
test_fit(EEG,'lsmr')
try
test_fit(EEG,'par-lsmr')
catch
    warning('parallel code had an error, no parallel toolbox?')
end
test_fit(EEG,'pinv')
test_fit(EEG,'matlab')
% test_fit(EEG,'glmnet') % the mex was crashing at the time

%% add some outliers that could throw off LSMR
% lsmr calculates tolerances based on the data. Thus even when data is not
% taken into account a huge outlier will completly mess up the fit because
% its not iteratively going deep enough. 
EEG.data(1,1) = 10^10;

test_fit(EEG,'lsmr')


function test_fit(EEG,method)
EEG = dc_glmfit(EEG,'method',method);
deviance = sum(sum(abs(squeeze(EEG.deconv.beta_dc(:,find(EEG.deconv.times>0,1),:)) - [3 2.5 -1.5; 3 2.5 -1.5])));
fprintf('deviance of %f for method %s \n',deviance,method)
assert(deviance<10^-5,'problem recovery beta values accurately enough')
