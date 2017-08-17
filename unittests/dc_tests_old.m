%%
display('This is a testing script, it will test SMOKE functionalities of the toolbox. I.e. it only detects if there will be a crash. no unit tests yet')
%%

init_unfold
%%
% set the seed so we get always the same result
oldseed = rng;
rng(0)

simCFG = [];
simCFG.nSignal  = 1:5;
simCFG.betalist = {3,0,0.5,[0 5 -2 3],[1,0.5]};
simCFG.betalist = {3,2,-1,[0 0 0 0],[0,0]};
simCFG.correlatedOverlap = {0 30 1 [0 0 0 0],0};
simCFG.data_length = 10;
EEG = simulate_data(simCFG);

simCFG.nSignal =1;
EEG2 = simulate_data(simCFG);
EEG.data = EEG.data + EEG2.data;

for e = 1:length(EEG.event)
    EEG.event(e).Y = 0;
end
for e = 1:length(EEG2.event)
    event = EEG2.event(e);
    event.latency = event.latency + randi(2000,1)-1000;
    if event.latency < 0
        event.latency = 0;
    end
    if event.latency > EEG.pnts
        event.latency = EEG.pnts;
    end
    event.type = 'stimulus';
    event.Y = 1;
    event.X_1 = 0;
    event.X_2 = 0;
    event.X_3 = 0;
    event.X_4 = 0;
    event.X_5 = 0;
    EEG.event(end+1) = event ;
    EEG
end

EEG = eeg_checkset(EEG,'eventconsistency');

ix = find(strcmp({EEG.event(:).type},'stimulus'));
for e = ix(10:15)
    EEG.event(e).Y = nan(1);
end

ix = find(strcmp({EEG.event(:).type},'trigger'));
for e = ix(10:15)
    EEG.event(e).X_1 = nan(1);
    EEG.event(e).X_2 = nan(1);
    EEG.event(e).X_3 = nan(1);
    EEG.event(e).X_4 = nan(1);
    EEG.event(e).X_5 = nan(1);
end
EEG = eeg_checkset(EEG,'eventconsistency');
display('simulation data generated')

rng(oldseed); %resets seed to user status
%%
display('Testing the generation of designmatrices')

for testCase = 1:7
    cfgDesign = [];
    cfgDesign.eventtype = {'trigger'};
    cfgDesign.coding = 'effects';
    switch testCase
        case 1
            cfgDesign.formula = 'y ~ 1'; % don't use '0' to drop intercept use '-1' % Wilkinson Notation
        case 2
            cfgDesign.formula = 'y ~ 1+ X_1+X_2+X_3'; % don't use '0' to drop intercept use '-1' % Wilkinson Notation
        case 3
            cfgDesign.formula = 'y ~ 1+ X_1+X_2+X_3'; % don't use '0' to drop intercept use '-1' % Wilkinson Notation
            cfgDesign.categorical = {'X_1','X_2'};

        case 4
            cfgDesign.formula = 'y ~ 1+ X_1+X_2+X_3+X_4'; % don't use '0' to drop intercept use '-1' % Wilkinson Notation
            cfgDesign.categorical = {'X_1','X_2','X_4'};
            cfgDesign.spline = {{'X_5',10}};
        case 5
            cfgDesign.formula = 'y ~ -1 + X_1+X_2+X_3+X_4'; % don't use '0' to drop intercept use '-1' % Wilkinson Notation
            cfgDesign.categorical = {'X_1','X_2','X_4'};
            cfgDesign.spline = {{'X_5',10}};
        case 6
            cfgDesign.eventtype = {{'trigger'},{'stimulus'}};
            cfgDesign.formula = {'y ~ 1 + X_1+X_2+X_3+X_4','y~1+Y'}; % don't use '0' to drop intercept use '-1' % Wilkinson Notation
            cfgDesign.categorical = {'X_1','X_2','X_4'};
            %cfgDesign.spline = {{'X_5',10}};
        case 7
            cfgDesign.eventtype = {{'trigger'},{'stimulus'}};
            cfgDesign.formula = {'y ~ 1 + X_1+X_2+X_3+X_4','y~1'}; % don't use '0' to drop intercept use '-1' % Wilkinson Notation
            cfgDesign.categorical = {'X_1','X_2','X_4'};
            %cfgDesign.spline = {{'X_5',10}};
            cfgDesign.coding = 'dummy';
    end


    EEG2 = dc_designmat(EEG,cfgDesign);
    fprintf('designmat test %i/%i without error (warnings are to be expected!) \n',testCase,7)
end
EEG = EEG2;





%%
%% Plotting tests

corr = dc_plotEventCorrmat(EEG,'eventtype',{'trigger'}); % correlation matrix
close(gcf)
dc_plotEventHistogram(EEG,'eventtype',{'trigger'})
close(gcf)

dc_plotDesignmat(EEG)
close(gcf)
display('basic plotting tests without error')

%%
%%
for testCase =  4%1:6
    cfgTimeshift = [];
    cfgTimeshift.timelimits = [-.1,0.8];
    switch testCase
        case 1
            cfgTimeshift.method = 'full';%'fourier';%'full';%'spxlines'
            %cfgTimeshift.timeshiftparam = 20; % either number of frequencies (param = * 2 because of cos/sin), or number of splines
            cfgTimeshift.sparse = 0;
        case 2
            cfgTimeshift.method = 'fourier';%'fourier';%'full';%'spxlines'
            cfgTimeshift.timeshiftparam = 20; % either number of frequencies (param = * 2 because of cos/sin), or number of splines
            cfgTimeshift.sparse = 0;
        case 3
            cfgTimeshift.method = 'splines';%'fourier';%'full';%'spxlines'
            cfgTimeshift.timeshiftparam = 20; % either number of frequencies (param = * 2 because of cos/sin), or number of splines
            cfgTimeshift.sparse = 0;
        case 4
            cfgTimeshift.method = 'full';%'fourier';%'full';%'spxlines'
            %cfgTimeshift.timeshiftparam = 20; % either number of frequencies (param = * 2 because of cos/sin), or number of splines
            cfgTimeshift.sparse = 1;
        case 5
            cfgTimeshift.method = 'fourier';%'fourier';%'full';%'spxlines'
            cfgTimeshift.timeshiftparam = 20; % either number of frequencies (param = * 2 because of cos/sin), or number of splines
            cfgTimeshift.sparse = 1;
        case 6
            cfgTimeshift.method = 'splines';%'fourier';%'full';%'spxlines'
            cfgTimeshift.timeshiftparam = 20; % either number of frequencies (param = * 2 because of cos/sin), or number of splines
            cfgTimeshift.sparse = 1;
    end
    EEG2 = dc_timeexpandDesignmat(EEG,cfgTimeshift);
    fprintf('timeshift test %i/%i without error \n',testCase,6)
end

EEG = EEG2;


%% Cleaning Test


%% fit the timeshifted DesignMat
EEG2= dc_glmfit(EEG,'channel',[1],'method','lsmr');
EEG3= dc_glmfit(EEG,'channel',[1],'method','pinv');
EEG4= dc_glmfit(EEG,'channel',[1],'method','matlab');
%EEG= dc_glmfit(EEG,'channel',[1],'method','par-lsmr');

reldiff = @(A,B)abs(A-B)/((max(A(:))+max(B(:)))/2);

eps = 10^-3; %tolerance of (0.1%)
diff = [reldiff(EEG2.deconv.dcBeta(:),EEG3.deconv.dcBeta(:)),... 
reldiff(EEG2.deconv.dcBeta(:),EEG4.deconv.dcBeta(:)),...
reldiff(EEG3.deconv.dcBeta(:),EEG4.deconv.dcBeta(:))];
fprintf('maximal difference between algorithms: %.3f%%\n',max(diff(:))*100)
if max(diff(:))>eps
    error('we have too large difference between algorithms!')
end


%%
