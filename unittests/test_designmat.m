function test_designmat()
% Multi-Test

testCase = [15]

EEGsim = simulate_test_case(testCase,'noise',0,'basis','box');
cfgDesign = [];

cfgDesign.coding = 'dummy';
cfgDesign.formula   = {'y~1',       'y~1+cat(conditionA)*continuousA', 'y~1+spl(splineA,5)+spl(splineB,5)+continuousA'};
cfgDesign.eventtypes = {'stimulus1', 'stimulus2',                       'stimulus3'};




uf_designmat(EEGsim,cfgDesign);


%%
% test some bugs that crashed the function (e.g. #4)
cfgDesign.formula{1} = 'xyz ~1';
uf_designmat(EEGsim,cfgDesign);
cfgDesign.formula{1} = 'a ~1';
uf_designmat(EEGsim,cfgDesign);

cfgDesign.formula{3} = 'a ~1+  spl(splineA ,5)'; % #3
uf_designmat(EEGsim,cfgDesign);
%% test two events with same predictorName
for e = 1:length(EEGsim.event)
    EEGsim.event(e).evt = e;
end
EEG2 = uf_designmat(EEGsim,'eventtypes',{'stimulus1' 'stimulus2'},'formula',{'y~evt','y~evt'});
assert(size(EEG2.unfold.X,2) == 4);


%% Only Interaction bug check
% We had a bug where specifying only the interaction without main effects
% jumbles things
EEG2 = uf_designmat(EEGsim,'eventtypes','stimulus2','formula','y~conditionA:continuousA');


assert(strcmp(EEG2.unfold.variablenames{EEG2.unfold.cols2variablenames(end)},'conditionA:continuousA'))

%% check higher order interactions
EEG2 = EEGsim;
for e =1:length(EEG2.event)
    EEG2.event(e).conditionB = randi(2,1);
end
EEG2.event = rmfield(EEG2.event,{'splineA','splineB'});
EEG2 = uf_designmat(EEG2,'eventtypes','stimulus2','formula','y~conditionB*conditionA:continuousA');

assert(size(EEG2.unfold.X,2) == 4)
assert(strcmp(EEG2.unfold.variablenames{EEG2.unfold.cols2variablenames(end)},'conditionA:conditionB:continuousA'))
%% Renaming checks
cfgDesign = [];
cfgDesign.coding = 'dummy';
cfgDesign.formula   = {'y~1+spl(splineA,4)+conditionA',       'y~1+conditionA*continuousA+splineA', 'y~1+spl(splineA,5)+spl(splineB,5)+continuousA'};
cfgDesign.eventtypes = {'stimulus1', 'stimulus2',                       'stimulus3'};

% fill in spline A for all events
for e= 1:length(EEGsim.event)
    
    EEGsim.event(e).splineA = rand(1);
end
EEGtmp = uf_designmat(EEGsim,cfgDesign);

shouldBeFunction(EEGtmp,{'(Intercept)','conditionA','splineA','2_(Intercept)','2_conditionA','continuousA','2_splineA','2_conditionA:continuousA','3_(Intercept)','3_continuousA','3_splineA','splineB'},'variablenames');


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
EEGtmp = uf_designmat(EEGsim,cfgDesign);



shouldBeFunction(EEGtmp,{'(Intercept)','conditionA','conditionB','2_(Intercept)','conditionC','2_conditionA:2_conditionB:conditionC','3_(Intercept)','3_conditionB','3_conditionC','3_conditionB:3_conditionC'},'variablenames')

%% Test the categorical rereferencing
EEGtest = EEGsim;
EEGtest.event = EEGtest.event(1:3:end);
catAlist = {'A','B','C'};
catBlist = {1,2,3};
for e= 1:length(EEGtest.event)
    
    EEGtest.event(e).catA= catAlist{mod(e,3)+1};
    EEGtest.event(e).catB= catBlist{mod(e,3)+1};
    
end

cfgDesign = [];
cfgDesign.coding = 'dummy';
cfgDesign.formula   = {'y~1+ cat(catA)+ catB'};
cfgDesign.eventtypes = {'stimulus1'};
cfgDesign.categorical = {'catA',{'B','C','A'};
                         'catB',{2,3,1}};
EEGtmp = uf_designmat(EEGtest,cfgDesign);


shouldBeFunction(EEGtmp,{'(Intercept)'  'catA_C'  'catA_A'  'catB_3'  'catB_1'},'colnames')
assert(all(strcmp('categorical',EEGtmp.unfold.variabletypes) == [0 1 1]))
%% Specifying only one reference category
cfgDesign.categorical = {'catB',{2,3,1}};
EEGtmp = uf_designmat(EEGtest,cfgDesign);
shouldBeFunction(EEGtmp,{'(Intercept)'  'catA_B'  'catA_C'  'catB_3'  'catB_1'},'colnames')
%% add a non-categorical variable
cfgDesign = [];
cfgDesign.eventtypes = {'stimulus1'};
cfgDesign.formula   = {'y~1+ catB + cat(catA)'};
EEGtmp = uf_designmat(EEGtest,cfgDesign);

shouldBeFunction(EEGtmp,{'(Intercept)','catA_B'  'catA_C','catB'},'colnames')

%% check that they are automatically sorted if not specified differently
cfgDesign = [];
cfgDesign.formula   = {'y~1+ cat(catA)+ cat(catB)'};
cfgDesign.eventtypes = {'stimulus1'};

EEGtest.event(1).catA = 'C';
EEGtest.event(1).catB = 3;
EEGtmp = uf_designmat(EEGtest,cfgDesign);

shouldBeFunction(EEGtmp,{'(Intercept)'  'catA_B'  'catA_C'  'catB_2'  'catB_3'},'colnames')

end

function shouldBeFunction(EEG,shouldBe,field)
for k = 1:length(EEG.unfold.(field))
    is = EEG.unfold.(field){k};
    assert(strcmp(is,shouldBe{k}),sprintf('error in %s, should be %s',is,shouldBe{k}))
end
end