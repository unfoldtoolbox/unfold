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


%% Only Interaction bug check
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
%% Renaming checks
cfgDesign = [];
cfgDesign.coding = 'dummy';
cfgDesign.formula   = {'y~1+spl(splineA,4)+conditionA',       'y~1+conditionA*continuousA+splineA', 'y~1+spl(splineA,5)+spl(splineB,5)+continuousA'};
cfgDesign.eventtypes = {'stimulus1', 'stimulus2',                       'stimulus3'};

% fill in spline A for all events
for e= 1:length(EEGsim.event)
    
EEGsim.event(e).splineA = rand(1);
end
EEGtmp = dc_designmat(EEGsim,cfgDesign);

shouldBe = {'(Intercept)','conditionA','splineA','2_(Intercept)','2_conditionA','continuousA','2_splineA','2_conditionA:continuousA','3_(Intercept)','3_continuousA','3_splineA','splineB'};
for k = 1:length(EEGtmp.deconv.variablenames)
    is = EEGtmp.deconv.variablenames{k};
    assert(strcmp(is,shouldBe{k}),sprintf('error in %s, should be %s',is,shouldBe{k}))
end

%% Renaming check interactions
cfgDesign = [];
cfgDesign.coding = 'dummy';
cfgDesign.formula   = {'y~1+ conditionA + conditionB','y~1+conditionC + conditionA:conditionB:conditionC','y~1+conditionC*conditionB'};
cfgDesign.eventtypes = {'stimulus1', 'stimulus2',                       'stimulus3'};

% fill in spline A for all events
for e= 1:length(EEGsim.event)
    
EEGsim.event(e).conditionA = randi(3)-1;
EEGsim.event(e).conditionB = randi(3)-1;
EEGsim.event(e).conditionC = randi(3)-1;
end
EEGtmp = dc_designmat(EEGsim,cfgDesign);


shouldBe = {'(Intercept)','conditionA','conditionB','2_(Intercept)','conditionC','2_conditionA:2_conditionB:conditionC','3_(Intercept)','3_conditionB','3_conditionC','3_conditionB:3_conditionC'};
for k = 1:length(EEGtmp.deconv.variablenames)
    is = EEGtmp.deconv.variablenames{k};
    assert(strcmp(is,shouldBe{k}),sprintf('error in %s, should be %s',is,shouldBe{k}))
end

end
