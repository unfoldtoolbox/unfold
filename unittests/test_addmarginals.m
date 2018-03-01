function [] = test_addmarginals()

%%

testCase = [16]
    
EEG = simulate_test_case(testCase,'noise',0,'basis','box');
cfgDesign = [];

cfgDesign.coding = 'dummy';
cfgDesign.formula   = 'y~1 + cat(conditionA) + continuousA'% + contB'% + spl(splineA,5)';
cfgDesign.eventtypes = 'stimulusA';

EEG = uf_designmat(EEG,cfgDesign);
EEG = uf_timeexpandDesignmat(EEG,'timelimits',[-0.5 1.5]);

EEG = uf_glmfit(EEG);

ufresult = uf_condense(EEG);

unfold3 = uf_predictContinuous(ufresult,'method','average');
ufresult = uf_predictContinuous(ufresult,'predictAt',{{'continuousA',[0 1]},{'contB',[0 1]},{'splineA',[-1 0 1]}});
unfold2 = uf_addmarginal(ufresult);



end

