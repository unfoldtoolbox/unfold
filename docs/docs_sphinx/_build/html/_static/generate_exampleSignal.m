if ~exist('f','var') || ~isvalid(f)
    f = figure
else
    clf(f)
end
nSubplots = 4;
fs = 300;
%%
ev1 = [0 4 10 14.5];
rng(4)
ev2 = cumsum(gamrnd(35,35,15,1)/1000)';
% ev2 = rand(15,1)*20+0.5;

subplot(nSubplots,1,1)
vline(ev1,'r')
vline(ev2,'g')

time = linspace(0,4,4*fs);
dipPos = hanning(600);
dipNeg = hanning(300);
resp1 = zeros(length(time),1);
resp1(1:(length(dipPos))) = resp1(1:length(dipPos)) + 1*dipPos;
resp1(500:(500+length(dipNeg)-1)) = resp1(500:(500+length(dipNeg)-1)) + -0.5 * dipNeg;

resp2 = zeros(length(time),1);
resp2(1:length(dipNeg)) = resp2(1:length(dipNeg)) + 0.5*dipNeg;
xlim([0,20]);set(gca,'box','off')
%%
subplot(nSubplots,1,2)
plot(time,resp1), hold all
plot(time,resp2)
xlim([0,20])
set(gca,'box','off')
%%
subplot(nSubplots,1,3)
hold all
for e = ev1
   plot(time+e,resp1,'r')
end

for e = ev2
   plot(time+e,resp2,'g')
end
xlim([0,20]);set(gca,'box','off')

%%
subplot(nSubplots,1,4)
t = 0:1/fs:(30-1/fs);
y = zeros(size(t))';
for e = ev1
    ix = get_min(e,t);
    y(ix:(ix+length(resp1)-1)) = y(ix:(ix+length(resp1)-1)) + resp1;
end

for e = ev2
    ix = get_min(e,t);
    y(ix:(ix+length(resp2)-1)) = y(ix:(ix+length(resp2)-1)) + resp2;
end
plot(t,y)
xlim([0,20]);set(gca,'box','off')