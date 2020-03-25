function [] = test_addmarginals()

%%

testCase = [15]
    
EEG = simulate_test_case(testCase,'noise',0,'basis','box');
cfgDesign = [];

cfgDesign.codingschema = 'reference';
cfgDesign.formula   = {'y~1 + cat(conditionA) + continuousA','y~1'};
cfgDesign.eventtypes = {{'stimulus2'},{'stimulus1','stimulus3'}};

EEG = uf_designmat(EEG,cfgDesign);
EEG = uf_timeexpandDesignmat(EEG,'timelimits',[-0.5 1.5]);

EEG = uf_glmfit(EEG);

ufresult = uf_condense(EEG);

ufresult_average = uf_predictContinuous(ufresult,'auto_method','average');
ufresult_predict = uf_predictContinuous(ufresult,'predictAt',{{'continuousA',[0 1]},{'contB',[0 1]},{'splineA',[-1 0 1]}});
ufresult_marginal = uf_addmarginal(ufresult_predict);



% the stimulus3 sould be unaffected, it's just an intercept
assert(1==near(ufresult_marginal.beta(:,:,5),ufresult_predict.beta(:,:,5))) % this column should not have been impacted.

% the stimulus2 conditionA should be sum of cont av + interaction +
% categorical predictor
assert(1==near(ufresult_marginal.beta(:,:,2),sum(ufresult_average.beta(:,:,[1 2 3]),3)))

% the stimulus 2, parameter 3 is continuous = 0, so it should be the same
% as the intercept
assert(1==near(ufresult_marginal.beta(:,:,3), ufresult.beta(:,:,1)))

% the stimulus 2, parameter 1 should be at continuous = mean(continuous)
% {that's 50} so it should be continuous-average + intercept
assert(1==near(ufresult_marginal.beta(:,:,1), sum(ufresult_average.beta(:,:,[1,3]),3)))


% Testing Bug #82
%%
EEG = simulate_test_case(15,'noise',0,'basis','box');
for e = 1:length(EEG.event)
    EEG.event(e).type= ['test_' num2str(randi(6,1))];
end
cfgDesign = [];
cfgDesign.codingschema = 'reference';
cfgDesign.formula   = {'y~1+continuousA','y~1+continuousA'};
cfgDesign.eventtypes = {{'test_1','test_5','test_3'},{'test_3','test_5'}};

EEG = uf_designmat(EEG,cfgDesign);
EEG = uf_timeexpandDesignmat(EEG,'timelimits',[-0.5 1.5]);

EEG = uf_glmfit(EEG);

ufresult = uf_condense(EEG);

ufresult_predict = uf_predictContinuous(ufresult,'predictAt',{{'continuousA',[0 1]},{'contB',[0 1]},{'splineA',[-1 0 1]}});

% this used to throw an error due to the event matching (because we want to
% add the marginal only to the event it belongs to).
ufresult_marginal = uf_addmarginal(ufresult_predict);
%% Testing AME and EME
rng(1)
    % Basis functions
intercept = struct();
intercept.eventname = 'stimulusA';
intercept.type = 'intercept';
intercept.predictorName = 'intercept';
intercept.overlap = 0; % note that 0, means that the AVERAGE overlap is exactly one signal-length (i.e. 1s in our simulations). That means if we do not want ANY overlap, we specify -1 here!
intercept.effectsize = 0;

cont = struct();
cont.eventname = 'stimulusA';
cont.type = 'continuous';
cont.overlap = 0;
cont.predictorName = 'continuousA';
cont.effectsize = 1;
cont.range = [0,100];

spline = struct();
spline.eventname = 'stimulusA';
spline.type = 'spline';
spline.overlap  = 0;
spline.predictorName = 'splineA';
spline.effectsize = nan;
spline.range = [-2,2];
spline.function = @(x)x.^3;

signals{1} = intercept;
signals{1}.range = nan;
signals{1}.function = nan;
signals{1}.overlap = -1;

signals{1}(2) = spline;
signals{1}(2).overlap = nan; % not yet defined
signals{1}(2).effectsize= 1;
signals{1}(2).predictorName= 'splineA';
signals{1}(2).function = @(x)(x+0.5).^3;
signals{1}(2).range = [-2,2];

signals{1}(3) = spline;
signals{1}(3).overlap = nan; % not yet defined
signals{1}(3).effectsize= 0.00;
signals{1}(3).predictorName= 'splineB';
signals{1}(3).range = [-2,2];
signals{1}(3).function = @(x)x.^3;

EEG = simulate_data(signals,'noise',0,'basis','box','datalength',500*60);
%%

cfgDesign = [];
cfgDesign.codingschema = 'reference';
cfgDesign.formula   ='y~1 +spl(splineA,5)+spl(splineB,5)';
cfgDesign.eventtypes = 'stimulusA';

EEG = uf_designmat(EEG,cfgDesign);
EEG = uf_timeexpandDesignmat(EEG,'timelimits',[-0.5 1.5]);

EEG = uf_glmfit(EEG);

ufresult = uf_condense(EEG);
%%
uf_plotParam((uf_predictContinuous(ufresult,'predictAt',{{'splineA',linspace(-2,2,11)}})))
%%
uf_plotParam(uf_addmarginal(uf_predictContinuous(ufresult,'predictAt',{{'splineA',linspace(-2,2,11)}}),'marginal','MEM'))
pAt = {{'splineA',linspace(-2,2,11)},{'splineB',linspace(-2,2,11)}};
uf_plotParam(uf_addmarginal(uf_predictContinuous(ufresult,'predictAt',pAt),'marginal','MEM'))
uf_plotParam(uf_addmarginal(uf_predictContinuous(ufresult,'predictAt',pAt),'marginal','AME'))
%%
x = uf_predictContinuous(ufresult,'predictAt',{{'splineA',0},{'splineB',0}});
assert(near(sum(x.beta(:,6,1:2)),0.5^3)==1)
%%
x = uf_addmarginal(uf_predictContinuous(ufresult),'marginal','MEM');
assert(near((mean([EEG.event.splineA])+0.5)^3,x.beta(:,6,1))==1);
%%
x = uf_addmarginal(uf_predictContinuous(ufresult),'marginal','AME');
%theoretical:
%mean((0.5+linspace(-2,2,100)).^3)
%practical:
assert(near(mean((([EEG.event.splineA])+0.5).^3),x.beta(:,6,1))==1);

end

