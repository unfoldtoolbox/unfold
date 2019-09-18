
%% Simulate
% Basis functions
intercept = struct();
intercept.eventname = 'stimulusA';
intercept.type = 'intercept';
intercept.predictorName = 'intercept';
intercept.overlap = 0; % note that 0, means that the AVERAGE overlap is exactly one signal-length (i.e. 1s in our simulations). That means if we do not want ANY overlap, we specify -1 here!
intercept.effectsize = 3;

singleFactor = struct();
singleFactor.eventname = 'stimulusA';
singleFactor.type = '1x2';
singleFactor.overlap = 0;
singleFactor.predictorName = 'conditionA';
singleFactor.effectsize = 1;

intercept = orderfields(intercept);
singleFactor= orderfields(singleFactor);


signals = [];
signals{1} = struct();
signals{1} = intercept;
signals{1}(1).overlap = 0.85;
signals{1}(2) = singleFactor;
signals{1}(2).overlap = 0.95;
signals{1}(2).effectsize= 0;
%%
cfg.plot = 0;
cfg.filterchoices = [0.01,0.1,1,10];
result.betas = [];
result.filter = [];
result.repeat = [];
for k = 1:20

EEG = simulate_data(signals,'basis','hanning','srate',100,'noise',0,'datalength',3000);

EEG.data = EEG.data(1,:)+(step(dsp.ColoredNoise(2,EEG.pnts)) +0.5*step(dsp.ColoredNoise(1,EEG.pnts)))';
EEG_raw = EEG;
%% Analysis
if cfg.plot
figure
subplot(3,1,1)
plot(EEG.times,EEG.data), xlabel('time [s]')
end

for filt = cfg.filterchoices
cfgDesign = [];
cfgDesign.eventtypes = {'stimulusA'};
cfgDesign.codingschema = 'effects';
cfgDesign.formula = 'y ~ 1 + conditionA';
cfgDesign.categorical = {'conditionA'};

EEG = pop_eegfiltnew(EEG_raw,filt,[]);
EEG = uf_designmat(EEG,cfgDesign);
% error
timelimits = [-0.8,1];
EEG = uf_timeexpandDesignmat(EEG,'timelimits',timelimits,'method','full','timeexpandparam',7);
EEG= uf_glmfit(EEG,'method','lsmr');

%
ufresult = uf_condense(EEG);
if cfg.plot
subplot(3,1,2)
hold all
plot(ufresult.times,ufresult.beta(1,:,1))
subplot(3,1,3)
hold all
plot(ufresult.times,ufresult.beta(1,:,2))
end

result.betas(end+1,:,:) = squeeze(ufresult.beta(1,:,:));
result.filter(end+1) = filt;
result.repeat(end+1) = k;

end
if cfg.plot
tmp = num2cell(cfg.filterchoices);
tmp = cellfun(@(x)num2str(x),tmp,'UniformOutput',0);
legend(tmp)
end
end

%%
g = gramm('x',ufresult.times,'y',result.betas(:,:,1),'color',result.filter);
% g.geom_line();
g.stat_summary()
g.draw();
%%
g = gramm('x',ufresult.times,'y',result.betas(:,:,2),'color',result.filter);
% g.geom_line()
g.stat_summary();
g.draw();