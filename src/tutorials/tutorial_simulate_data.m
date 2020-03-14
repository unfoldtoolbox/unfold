function EEG = tutorial_simulate_data(type)

switch type
    case '2x2'
        % Basis functions
        sprev = rng(3); %set random-generator seed
        
        signals{1}= struct();
        signals{1}.eventname = 'fixation';
        signals{1}.type = 'intercept';
        signals{1}.predictorName = 'intercept';
        signals{1}.overlap = 0.5; % note that 0, means that the AVERAGE overlap is exactly one signal-length (i.e. 1s in our simulations). That means if we do not want ANY overlap, we specify -1 here!
        signals{1}.effectsize = 3;
        
        
        signals{1}(2).eventname = 'fixation';
        signals{1}(2).type = '1x2';
        signals{1}(2).overlap = -0.5;
        signals{1}(2).predictorName = 'stimulusType';
        signals{1}(2).effectsize = 1.5;
        
        signals{1}(3).eventname = 'fixation';
        signals{1}(3).type = '1x2';
        signals{1}(3).overlap = 0;
        signals{1}(3).predictorName = 'color';
        signals{1}(3).effectsize = -0.5;
        
        EEG = simulate_data(signals,'srate',200,'basis','hanning');
        
        rng(sprev) % reset random generator to previous seed
        EEG.times = []; 
        EEG = eeg_checkset(EEG);
        % renaming events for better interpretability
        for e = 1:length(EEG.event)
            if EEG.event(e).stimulusType == 1
                EEG.event(e).stimulusType = 'car';
            else
                EEG.event(e).stimulusType= 'face';
            end
            if EEG.event(e).color == 1
                EEG.event(e).color= 'green';
            else
                EEG.event(e).color= 'red';
            end
        end
    case 'trf'
        % makes use of the 2x2 and adds a predictor
        EEG = tutorial_simulate_data('2x2');
        % Add a "channels" reflecting our TRF predictor values
        EEG.data(end+1,:) = rand(size(EEG.data,2),1);
        
        % generate a kernel that we want to recover
        kernel = hanning(100);
        kernel(50:end) = kernel(50:end) - 0.7*hanning(51);
        tmp  = conv(EEG.data(end,:),kernel);
        tmp = tmp - mean(tmp);
        
        % add the TRF to the already simulated EEG data
        EEG.data(1,:) = EEG.data(1,:) + tmp(1:size(EEG.data,2));
    case 'erpimage'
         sprev = rng(3); %set random-generator seed
        
        signals{1}= struct();
        signals{1}.eventname = 'stimulus';
        signals{1}.type = 'intercept';
        signals{1}.predictorName = 'intercept';
        signals{1}.overlap = -1; % note that 0, means that the AVERAGE overlap is exactly one signal-length (i.e. 1s in our simulations). That means if we do not want ANY overlap, we specify -1 here!
        signals{1}.effectsize = 2;
        
        
        signals{1}(2).eventname = 'stimulus';
        signals{1}(2).type = '1x2';
        signals{1}(2).overlap = -0.5;
        signals{1}(2).predictorName = 'stimulusType';
        signals{1}(2).effectsize = 1.5;
        
        signals{2}(1).eventname = 'buttonpress';
        signals{2}(1).type = 'intercept';
        signals{2}(1).overlap = 0.5;
        signals{2}.predictorName = 'intercept';
        signals{2}(1).effectsize = 2.5;
        
        EEG = simulate_data(signals,'srate',200,'basis','posneg');
        
        rng(sprev) % reset random generator to previous seed
        EEG.times = []; 
        EEG = eeg_checkset(EEG);
        % renaming events for better interpretability
        for e = 1:length(EEG.event)
            if EEG.event(e).type == "stimulus"
            if EEG.event(e).stimulusType == 1
                EEG.event(e).stimulusType = 'car';
            else
                EEG.event(e).stimulusType= 'face';
            end
            end
        end
        
end