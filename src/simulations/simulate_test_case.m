function EEG = simulate_test_case(sim,varargin)

sprev = rng(1);
assert(ismember(sim,1:15),'requested signal has to be between 1-14')
% Varargin is loopt to simulate_data2

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



%%

% OVERLAP of 1 means NO OVERLAP for the condition effects
signals = [];
signals{1} = struct();
switch sim
    
    case 1 % One Intercept only no overlap
        
        signals{1} = intercept;
        signals{1}(1).overlap = -1;
        
    case 2 % One Intercept only with overlap
        signals{1} = intercept;
        signals{1}(1).overlap = 0.4;
        
    case 3 % One 1x2 Factor with an effect but not on timing
        signals{1} = intercept;
        signals{1}(1).overlap = -1;
        
        signals{1}(2) = singleFactor;
        signals{1}(2).effectsize= 1;
    case 4 % One 1x2 Factor with no effect except on timing
        signals{1} = intercept;
        signals{1}(1).overlap = 0.85;
        signals{1}(2) = singleFactor;
        signals{1}(2).overlap = 0.95;
        signals{1}(2).effectsize= 0;
        
    case 5 % 	- One continuous Factor with abs(slope) > 0 but not on timing
        signals{1} = intercept;
        signals{1}.range = nan;
        signals{1}.overlap = -1;
        signals{1}(2) = cont;
        signals{1}(2).overlap = 0;
        signals{1}(2).effectsize= 0.02;
    case 6 % 	- One continuous Factor with 0-slope but slope on timing
        signals{1} = intercept;
        signals{1}.range = nan;
        signals{1}(2) = cont;
        signals{1}(2).overlap = 0.4;
        signals{1}(2).effectsize= 0;
    case 7 % intercept plus simple cubic effect, no overlap
        signals{1} = intercept;
        signals{1}.range = nan;
        signals{1}.function = nan;
        signals{1}.overlap = -1;
        
        signals{1}(2) = spline;
        signals{1}(2).overlap = nan; % not yet defined
        signals{1}(2).effectsize= 0.01;
    case 8 % intercept plus simple cubic effect, with overlap
        signals{1} = intercept;
        signals{1}.range = nan;
        signals{1}.function = nan;
        signals{1}.overlap = 1;
        
        signals{1}(2) = spline;
        signals{1}(2).overlap = nan; % not yet defined
        signals{1}(2).effectsize= 0.01;
        
    case 9 % 	- Two events, no overlap
        signals{1} = intercept;
        signals{1}.overlap = -1;
        signals{2} = intercept;
        signals{2}.overlap = -1;
        signals{2}.eventname = 'stimulusB';
    case 10 % 	- Two events, with overlap
        signals{1} = intercept;
        signals{1}.overlap = 0.4;
        signals{2} = intercept;
        signals{2}.overlap = 0.4;
        signals{2}.eventname = 'stimulusB';
        
    case 11 % 	- Three events - Simple Intercepts, no overlap
        signals{1} = intercept;
        signals{1}.overlap = -1;
        signals{2} = intercept;
        signals{2}.overlap = -1;
        signals{2}.eventname = 'stimulusB';
        signals{3} = intercept;
        signals{3}.overlap = -1;
        signals{3}.eventname = 'stimulusC';
    case 12 % 	- Three events - Simple Intercepts, with overlap
        signals{1} = intercept;
        signals{1}.overlap = 0.4;
        signals{2} = intercept;
        signals{2}.overlap = 0.4;
        signals{2}.eventname = 'stimulusB';
        signals{3} = intercept;
        signals{3}.overlap = 0.4;
        signals{3}.eventname = 'stimulusC';
    case 13 %   - Three events - Different formulas for each, intercept only, 1x2 anova, continuous predictor, no overlap
        signals{1} = intercept;
        signals{1}.overlap = -1;
        
        signals{2} = intercept;
        signals{2}.eventname = 'stimulusB';
        signals{2}(1).overlap = -1;
        signals{2}(2) = singleFactor;
        signals{2}(2).effectsize= 1;
        
        
        signals{3} = intercept;
        signals{3}.eventname = 'stimulusC';
        signals{3}(1).overlap = -1;
        signals{3}(1).range = nan;
        
        signals{3}(2) = cont;
        signals{3}(2).overlap = 0;
        signals{3}(2).effectsize= 0.02;
    case 14 %   - Three events - Different formulas for each, intercept only, 1x2 anova, continuous predictor, with overlap
        signals{1} = intercept;
        signals{1}.overlap = 0.4;
        signals{2} = intercept;
        signals{2}(1).eventname = 'stimulusB';
        signals{2}(2) = singleFactor;
        signals{2}(2).overlap = 0.8;
        signals{2}(2).effectsize= 1;
        
        
        signals{3} = intercept;
        signals{3}.eventname = 'stimulusC';
        signals{3}(1).range = nan;
        signals{3}(2) = cont;
        signals{3}(2).overlap = 0.8;
        signals{3}(2).effectsize= 0.02;
        
    case 15 %  The overkill ;)
        signals{1} = intercept;
        signals{1}.overlap = 0.4;
        signals{2} = intercept;
        
        signals{2}(1).eventname = 'stimulusB';
        signals{2}(2) = singleFactor;
        signals{2}(1).range = nan;

        signals{2}(2).overlap = 0.8;
        signals{2}(2).effectsize= 1;
        signals{2}(2).range = nan;
        signals{2}(3) = cont;
        signals{2}(3).overlap = 0.8;
        signals{2}(3).effectsize= 0.02;
        
        
        signals{3} = intercept;
        signals{3}.range = nan;
        signals{3}.function = nan;
        signals{3}.overlap = -1;
        
        signals{3}(2) = spline;
        signals{3}(2).overlap = nan; % not yet defined
        signals{3}(2).effectsize= 0.01;
        signals{3}(2).predictorName = 'splineA';
        signals{3}(3) = spline;
        signals{3}(3).overlap = nan; % not yet defined
        signals{3}(3).effectsize= 0.01;
        signals{3}(3).predictorName = 'splineB';
        
        cont2 = cont;
        cont2.function = nan;
        signals{3}(4) = cont2;
        signals{3}(4).overlap = 0.2; % not yet defined
        signals{3}(4).effectsize= 0.02;
        
        signals{1}(1).eventname = 'stimulus1';
        signals{2}(1).eventname = 'stimulus2';
        signals{3}(1).eventname = 'stimulus3';
        
    otherwise
        error('unknown simulation')
end
EEG = simulate_data(signals,varargin{:});
rng(sprev);