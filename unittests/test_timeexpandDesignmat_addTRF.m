function test_uf_timeexpandedDesignmat_addTRF()
%% 
% Generates Data in EEG.data(1,:) which is the time-shifted sum of two
% convolutions of two random signals with a hanning kernel each.
% Recovers those kernels using TRF

EEGsim = eeg_emptyset();
EEGsim.data(2,1:1000) = rand(1000,1);
EEGsim.data(3,1:1000) = rand(1000,1);
kernel = hanning(15);
tmp = conv(EEGsim.data(2,:),kernel);
tmp2 = conv(EEGsim.data(3,:),-0.5*kernel);
EEGsim.data(1,:) = tmp(1:size(EEGsim.data,2));
EEGsim.data(1,:) = EEGsim.data(1,:) + tmp2(1:size(EEGsim.data,2));
EEGsim.srate = 1000;
EEGsim = eeg_checkset(EEGsim);

EEGsim2 = uf_timeexpandDesignmat_addTRF(EEGsim,'channel',2,'name','covA','timelimits',[-0.01,0.02]);
EEGsim2 = uf_timeexpandDesignmat_addTRF(EEGsim2,'channel',3,'name','covB','timelimits',[-0.01,0.02]);
EEGsim2 = uf_glmfit(EEGsim2);

assert(near(kernel,squeeze(EEGsim2.unfold.beta_dc(1,find(EEGsim2.unfold.times==0):find(EEGsim2.unfold.times==0)+14,1))')==1)
assert(near(-0.5*kernel,squeeze(EEGsim2.unfold.beta_dc(1,find(EEGsim2.unfold.times==0):find(EEGsim2.unfold.times==0)+14,2))')==1)