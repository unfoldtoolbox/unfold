function test_designmat()
% Multi-Test

testCase = [15]
    
EEGsim = simulate_test_case(testCase,'noise',0,'basis','box');
cfgDesign = [];

cfgDesign.coding = 'dummy';
cfgDesign.formula   = {'y~1',       'y~1+cat(conditionA)*continuousA', 'y~1+spl(splineA,5)+spl(splineB,5)+continuousA'};
cfgDesign.eventtypes = {'stimulus1', 'stimulus2',                       'stimulus3'};




dc_designmat(EEGsim,cfgDesign);


%%
% test some bugs that crashed the function (e.g. #4)
cfgDesign.formula{1} = 'xyz ~1';
dc_designmat(EEGsim,cfgDesign);
cfgDesign.formula{1} = 'a ~1';
dc_designmat(EEGsim,cfgDesign);

cfgDesign.formula{3} = 'a ~1+  spl(splineA ,5)'; % #3
dc_designmat(EEGsim,cfgDesign);
%% test two events with same predictorName
for e = 1:length(EEGsim.event)
   EEGsim.event(e).evt = e; 
end
EEG2 = dc_designmat(EEGsim,'eventtypes',{'stimulus1' 'stimulus2'},'formula',{'y~evt','y~evt'});
assert(size(EEG2.deconv.X,2) == 4);
%%
% We had a bug where specifying only the interaction without main effects
% jumbles things
EEG2 = dc_designmat(EEGsim,'eventtypes','stimulus2','formula','y~conditionA:continuousA');


assert(strcmp(EEG2.deconv.variablenames{EEG2.deconv.cols2variablenames(end)},'conditionA:continuousA'))

%% check higher order interactions
EEG2 = EEGsim;
for e =1:length(EEG2.event)
    EEG2.event(e).conditionB = randi(2,1);
end
EEG2.event = rmfield(EEG2.event,{'splineA','splineB'});
EEG2 = dc_designmat(EEG2,'eventtypes','stimulus2','formula','y~conditionB*conditionA:continuousA');

assert(size(EEG2.deconv.X,2) == 4)
assert(strcmp(EEG2.deconv.variablenames{EEG2.deconv.cols2variablenames(end)},'conditionA:conditionB:continuousA'))
%%
cfgDesign = [];
cfgDesign.coding = 'dummy';
cfgDesign.formula   = {'y~1+spl(splineA,4)+conditionA',       'y~1+cat(conditionA)*continuousA', 'y~1+spl(splineA,5)+spl(splineB,5)+continuousA'};
cfgDesign.eventtypes = {'stimulus1', 'stimulus2',                       'stimulus3'};

% fill in spline A for all events
for e= 1:length(EEGsim.event)
    
EEGsim.event(e).splineA = rand(1);
end
EEGtmp = dc_designmat(EEGsim,cfgDesign);
assert(strcmp(EEGtmp.deconv.variablenames{4},'2_(Intercept)'))
assert(strcmp(EEGtmp.deconv.variablenames{5},'2_conditionA'))
assert(strcmp(EEGtmp.deconv.variablenames{6},'continuousA'))
assert(strcmp(EEGtmp.deconv.variablenames{8},'3_(Intercept)'))
assert(strcmp(EEGtmp.deconv.variablenames{9},'3_continuousA'))
assert(strcmp(EEGtmp.deconv.variablenames{10},'3_splineA'))
assert(strcmp(EEGtmp.deconv.variablenames{11},'splineB'))
