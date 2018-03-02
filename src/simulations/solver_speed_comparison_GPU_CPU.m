% Load any arbitrary large EEG after the uf_timeexpand step
cfg = [];
cfg.sizemult = [1 4];
cfg.sizemult = [1 1];

X = EEG.unfold.Xdc;
X = repmat(X,cfg.sizemult);
y = double(EEG.data(1,:)');



tic
Xdcg = gpuArray(X);
yg = gpuArray(y);
tGPUinit = toc;


%% CPU VS GPU, LSQR vs LSMR
tic
cpu_blsqr = lsqr(X,y,10^-8,300);
tCPUsolve = toc;


tic
gpu_blsqr = lsqr(Xdcg,yg,10^-8,300);
tGPUsolve = toc;

tic
[cpu_blsmr,ISTOP,ITN] = lsmr(X,(y),[],10^-8,10^-8,[],300); % ISTOP = reason why algorithm has terminated, ITN = iterations
display(ITN)
tCPUlsmr = toc;



fprintf('tGPUinit:%.2f\n',tGPUinit)
fprintf('CPU-lsqr:%.2f\nGPU-lsqr:%.2f\nCPU-lsmr:%.2f\n',tCPUsolve,tGPUsolve,tCPUlsmr)
max(abs(gpu_blsqr - cpu_blsqr))
max(abs(cpu_blsqr - cpu_blsmr))

%CPU lsmr seems to be the fastest (15% faster than lsqr)
%GPU lsqr seems to be 2x slower than CPU lsmr


%% LSQR give the solution to the next step
% This is in general slightly slower than normal LSQR (not sure why though)
tList = [];
% solve one channel, give the solution to the next
for k = 1:30
first_solve = lsqr(X,double(EEG.data(k,:)'),10^-8,300,[],[]);
y = double(EEG.data(k+1,:)');
 
tic
first_solve = lsqr(X,y,10^-8,300,[],[],first_solve);
tSecondSolve= toc;

tic
first_solve = lsqr(X,y,10^-8,300);
tFirstSolve= toc;

fprintf('First:%.2fs \t second:%.2fs \t diff:%.2fs \n',tFirstSolve,tSecondSolve,tFirstSolve-tSecondSolve)

tList(end+1,1:2) = [tFirstSolve,tSecondSolve];
end

%% LSMR (hacked) give the solution to the next step
% This is in general slightly slower than normal LSQR (not sure why though)
tList = [];

% precondition
s= sqrt(sum(X.^2));
Xpre = X./s;

% solve one channel, give the solution to the next
for k = 1:30
[first_solve,ISTOP,ITN] = lsmr(Xpre,double(EEG.data(k,:)'),[],10^-8,10^-8,[],300); % ISTOP = reason why algorithm has terminated, ITN = iterations
y = double(EEG.data(k+1,:)');
 
tic
[b_pre,ISTOP,ITN] = lsmr(Xpre,(y),[],10^-8,10^-8,[],300,[],[],first_solve); % ISTOP = reason why algorithm has terminated, ITN = iterations

tSecondSolve= toc;

tic
[b_pre,ISTOP,ITN] = lsmr(Xpre,(y),[],10^-8,10^-8,[],300); % ISTOP = reason why algorithm has terminated, ITN = iterations
tFirstSolve= toc;

fprintf('First:%.2fs \t second:%.2fs \t diff:%.2fs \n',tFirstSolve,tSecondSolve,tFirstSolve-tSecondSolve)

tList(end+1,1:2) = [tFirstSolve,tSecondSolve];
end

%% Preconditioning on X
% about 2x as fast
s= sqrt(sum(X.^2));
Xpre = X./s;

tic
[b_pre,ISTOP,ITN] = lsmr(Xpre,(y),[],10^-8,10^-8,[],300); % ISTOP = reason why algorithm has terminated, ITN = iterations
b_pre = b_pre./s;
display(ITN)
tPre = toc;


tic
[b_nopre,ISTOP,ITN] = lsmr(X,(y),[],10^-8,10^-8,[],300); % ISTOP = reason why algorithm has terminated, ITN = iterations
display(ITN)
tnoPre = toc;

%%

tic 
cpu_blsqr = lsqr(Xpre,y,10^-8,300);
toc


tic 
cpu_blsqr = lsqr(X,y,10^-8,300);
toc

