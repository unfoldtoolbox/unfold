function test_timeexpandDesignmat()

testCase = [1]
    
EEGsim = simulate_test_case(testCase,'noise',0,'basis','box');
cfgDesign = [];

cfgDesign.coding = 'dummy';
cfgDesign.formula   = {'y~1'};
cfgDesign.eventtype = {'stimulusA'};




EEGsim = dc_designmat(EEGsim,cfgDesign);
%%
EEG2 = dc_timeexpandDesignmat(EEGsim,'timelimits',[-1 2])


events = sum(strcmp({EEG2.event(:).type},'stimulusA'));
timeshiftedvals = sum(EEG2.deconv.Xdc(:,1));
assert(events==timeshiftedvals)
% We had a bug where only negative entries did not timeexpand it
EEG2 = EEGsim;
EEG2.deconv.X = -EEG2.deconv.X;
EEG2 = dc_timeexpandDesignmat(EEG2,'timelimits',[-1 2])

events = sum(strcmp({EEG2.event(:).type},'stimulusA'));
timeshiftedvals = sum(EEG2.deconv.Xdc(:,1));
assert(events==-timeshiftedvals)

