function test_checkmodelfit
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
EEG.data(2,:,:) = EEG.data(1,:,:)+1000*(rand(size(EEG.data(1,:,:)))-.5); % add a second channel
%% Create the model
cfgDesign = [];
cfgDesign.codingschema = 'reference';
cfgDesign.formula   = {'y~1+continuousA+continuousB'};
cfgDesign.eventtypes = {'stimulusA'};

EEG = uf_designmat(EEG,cfgDesign);
EEG = eeg_checkset(EEG);
EEG = uf_timeexpandDesignmat(EEG,'timelimits',[0 1]);
EEG = uf_glmfit(EEG);

ufresult = uf_condense(EEG);
%%
r2 = uf_checkmodelfit(EEG,'method','R2','fold_event','pause','channel',1);

%%
r2 = uf_checkmodelfit(EEG,'method','R2','fold_event','pause');
cv_r2 = uf_checkmodelfit(EEG,'method','crossValR2','fold_event',{'pause'});
median(cv_r2,2)

%%
r2 = uf_checkmodelfit(EEG,'method','partialR2');
assert(size(r2,1) == 6)
r2_subset = uf_checkmodelfit(EEG,'method','partialR2','variablename',{'continuousA','continuousB'});
assert(size(r2_subset,1) == 4)

r2 = uf_checkmodelfit(EEG,'method','partialR2','channel',2);
assert(all(r2.channel == 2)) %108 channel names not coming through
%%
cv_r2 = uf_checkmodelfit(EEG,'method','crossValpartialR2','fold_event',{'pause'});
groupsummary(cv_r2,{'variablenames','channel'},{'median','std'},{'r2_ca'})
r2 = uf_checkmodelfit(EEG,'method','crossValpartialR2','variablename',{'continuousA','continuousB'},'fold_event','pause');
assert(length(unique(r2.variablenames))==2) % check variable selection
%% bug 108 Splines not working with partialR2
cfgDesign = [];
cfgDesign.codingschema = 'reference';
cfgDesign.formula   = {'y~1+spl(continuousA,10)+continuousB'};
cfgDesign.eventtypes = {'stimulusA'};

EEG = uf_designmat(EEG,cfgDesign);
EEG = uf_timeexpandDesignmat(EEG,'timelimits',[0 1]);
EEG = uf_glmfit(EEG);

ufresult = uf_condense(EEG);
%%
r2 = uf_checkmodelfit(EEG,'method','partialR2','fold_event','pause','channel',1);
assert(size(r2,1)==3)
r2 = uf_checkmodelfit(EEG,'method','crossValpartialR2','fold_event','pause');
assert(size(r2,1) == 2*3*10) % 2chan x 3 variables x 10 folds

