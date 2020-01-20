function test_glmfit
%%
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
signals{1}.overlap = 0.5;
signals{1}(2) = cont;
signals{1}(2).overlap = 0.5;
signals{1}(2).effectsize= 2.5;
signals{1}(3) = cont;
signals{1}(3).predictorName = 'continuousB';
signals{1}(3).overlap = 0.5;
signals{1}(3).effectsize= -1.5;

for k = 1:10
    EEG_tmp = simulate_data(signals,'noise',100,'datalength',20);
    EEG_tmp.data(1,end:end+1000) = 0;
    EEG_tmp.pnts = EEG_tmp.pnts+1000;
    if k == 1
        EEG = EEG_tmp;
    else
        EEG = pop_mergeset(EEG,EEG_tmp);
        
    end
    % add a "break" event
    
    EEG.event(end+1) = EEG.event(end);
    EEG.event(end).latency = EEG.pnts;
    EEG.event(end).type = 'pause';
    
end

%% Create the model
cfgDesign = [];
cfgDesign.coding = 'dummy';
cfgDesign.formula   = {'y~1+continuousA+continuousB'};
cfgDesign.eventtypes = {'stimulusA'};

EEG = uf_designmat(EEG,cfgDesign);
EEG = eeg_checkset(EEG);
EEG = uf_timeexpandDesignmat(EEG,'timelimits',[0 1]);
EEG = uf_glmfit(EEG);

EEG = uf_glmfit(EEG,'method','glmnet','fold_event',{'pause'});
ufresult = uf_condense(EEG);
%%
r2 = uf_modelcheck(EEG,'method','R2','fold_event','pause');
cv_r2 = uf_modelcheck(EEG,'method','crossValR2','fold_event',{'pause'});
r2_cv  = mean(cv_r2)
r2

%%
r2 = uf_modelcheck(EEG,'method','commonalityR2')

%%
r2 = uf_modelcheck(EEG,'method','crossValcommonalityR2','fold_event','pause');
groupsummary(r2,'variablenames',{'median','std'},{'r2_ca'})
