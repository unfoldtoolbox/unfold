function [EEG] = simulate_data(signals,varargin)

% 'eventname', 'stimulusA', 'type','intercept', 'overlap',0 [range 0 1],

% 'eventname', 'stimulusA', 'type','1x2', 'overlap',0 [range 0 1], 'eventname','factorA', 'effectsize', 1 (a basis of 1, +1 for if the effect is active)

% 'eventname', 'stimulusA', 'type','continuous', 'overlap',0 [correlation
% range 0 1], 'eventname','continuousA', 'effectsize', 1 (slope), 'range', [0,100]

% 'eventname', 'stimulusA', 'type','spline', 'overlap',0 [not sure how to do this], 'eventname','splineA', 'range', [-10,10], 'function', @(x)x^3

simCFG= finputcheck(varargin,...
    {'datalength','integer',[],10*60;
    'noise','real',[],1;
    'srate','integer',[],10;
    'basis','string',{'box','hanning','dirac','posneg'},'box'
    },'mode','ignore');

assert(~ischar(simCFG),simCFG)

% mean interstimulus interval
cfg = [];
cfg.srate = simCFG.srate;%Hz
cfg.pnts= simCFG.datalength*cfg.srate; % default: 10 times 60s, thus 10minute of data
cfg.noise = simCFG.noise;


%% ###############################################################
% create response function for 3 events

sig = struct();
sig.time= 0:1/cfg.srate:(1-1/cfg.srate); % 1 second stimulus
% sig.time = 0:1/cfg.srate:(0.2-1/cfg.srate);%
% warning('BEHINGER MODIFIED THIS')


sig.shape=zeros(1,length(sig.time),1);
switch simCFG.basis
    case 'box'
        sig.shape= ones(1,length(sig.time));
    case 'hanning'
        sig.shape = hanning(length(sig.time)); %P3
    case 'dirac'
        sig.shape = [1 zeros(1,length(sig.time)-1)];
    case 'posneg'
        sig.shape = [hanning(floor(length(sig.time)/2)); -hanning(ceil(length(sig.time)/2))]; %P3
end

% sig.shape= [sig.shape; zeros(1,80)'];


%% ##################

%% function to generate one row of X for a single event
    function X = generateX(signals)
        X = nan(1,length(signals));
        for variable = 1:length(signals)
            
            switch signals(variable).type
                case 'intercept'
                    X(variable) = 1;
                case '1x2'
                    X(variable) = binornd(1,0.5);
                case 'continuous'
                    r = signals(variable).range;
                    X(variable) = rand(1)*(r(2)-r(1))+r(1);
                case 'spline'
                    r = signals(variable).range;
                    X(variable) = rand(1)*(r(2)-r(1))+r(1);
                    %X(variable) = signals(variable).function(X(variable));
                otherwise
                    error('unknown variable type: %s', signals(variable).type)
            end
        end
        
    end

    function O = generateO(signals,oneX)
        O = 0;
        for variable = 1:length(signals)
            ov = (1-signals(variable).overlap);
            switch signals(variable).type
                case 'intercept'
                    O = O + 1000 * ov; %the basis
                case '1x2'
                    O = O + 1000 * ov*oneX(variable); % in addition to that (oneX has to be 0 or 1
                case 'continuous'
                    r = signals(variable).range;
                    O = O + 1000 *  ov * oneX(variable)/(r(2)-r(1));
                    
                case 'spline'
                    if ~isnan(ov)
                        warning('spline overlap-bias is not implemented yet')
                    end
                otherwise
                    error('unknown variable type: %s', signals(variable).type)
            end
        end
        
    end



timeToNextEvent = [];
eventSignal = [];
X = cell(size(signals));
runner = 1;
breakFlag = 0;
while sum(timeToNextEvent)<cfg.pnts
    if breakFlag == 1
        break
    end
    for sIx = 1:length(signals) % in case of multiple events
        XoneRow = generateX(signals{sIx});
        
        % Save the "designMatrix"
        if isempty(X{sIx})
            X{sIx} = XoneRow;
        else
            X{sIx}(end+1,:) = XoneRow;
        end
        % generate overlap (i.e. when the next signal should occur)
        % actually this generates the mean of the distribution of distance.
        % if this mean is very small, lots of overlap
        try
        XrowBefore = X{sIx}(end-1,:);
        catch
            XrowBefore = XoneRow;
        end
        overlap = generateO(signals{sIx},XrowBefore);
        
        m = overlap;
        
        if m <= 0
%            warning('mean of log-distributed distances was smaller equal 0, forcing to be at least 1ms: %f',m)
            m = 1;
        end
        
        v = 100^2;
        % taken from 'lognstat' matlab help:
        mu = log((m^2)/sqrt(v+m^2));
        sigma = sqrt(log(v/(m^2)+1));
        % to visualize adapt accordingly:
        % figure,hist([lognrnd(mu,sigma,1,1000)*0.01;lognrnd(9.5,0.5,1,1000)*0.01]',100)
        nextTrigger = lognrnd(mu,sigma,1,1)/1000*cfg.srate;
        
        t = ceil(nextTrigger);
        
        if (sum(timeToNextEvent)+t)>cfg.pnts
            breakFlag = 1;
            break
        end
        
        %fprintf('%i from %i \n',sum(timeToNextEvent),cfg.pnts)
        runner = runner+1;
        eventSignal = [eventSignal sIx];
        timeToNextEvent = [timeToNextEvent t];
    end
end
% this gives the event times, in "eventSignal" is written, which times
% belong to which events
eventTimes = cumsum(timeToNextEvent);

for sIx = 1:length(signals)
    
    % generate the simulated data
    EEG_tmp=generate_signal(X{sIx},eventTimes(eventSignal==sIx),sig,signals{sIx},cfg);
    
    % concatenate multiple signals
    if sIx == 1
        EEG = EEG_tmp;
        EEG.sim.signals = signals;
    else
        % Combine events (can have different fields which makes
        % concatenation difficult)
        fnB = fieldnames(EEG_tmp.event);
        fnA = fieldnames(EEG.event);
        
        fnAinB = find(~ismember(fnA,fnB));
        fnBinA = find(~ismember(fnB,fnA));
        
        ev1= EEG.event;
        ev2 = EEG_tmp.event;
        ev1 = rmfield(ev1,fnA(fnAinB));
        ev2 = rmfield(ev2,fnB(fnBinA));
        ev = [ev1;ev2];
        % Now fill in all of A, that was not in B
        for e =1:length(ev1)
            for f = fnAinB'
             ev(e).(fnA{f}) = EEG.event(e).(fnA{f});
            end
        end
        for e =1:length(ev2)
            for f = fnBinA'
             ev(e+length(EEG.event)).(fnB{f}) = EEG_tmp.event(e).(fnB{f});
            end
        end
        
        EEG.event = ev;%[EEG.event;EEG_tmp.event];
        % Combine Simulation
        EEG.sim.X = X;
        EEG.sim.separateSignal = [EEG.sim.separateSignal; EEG_tmp.sim.separateSignal];
        EEG.sim.noOverlapSignal = [EEG.sim.noOverlapSignal; EEG_tmp.sim.noOverlapSignal];
        EEG.eventTimes = eventTimes;
        % Combine Data
        EEG.data = EEG.data + EEG_tmp.data;
        EEG = eeg_checkset(EEG,'eventconsistency');
    end
end

end
%%
