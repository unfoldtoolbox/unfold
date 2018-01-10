rng(1); %speficy seed
b = [0.1 0.5 1 2 0.2 4 2]/2;
N = 200;

x = linspace(0.5,10.5,N)';
y = [];
yD = log(x) + randn(N,1)*0.1;
yD = yD - min(yD);
% Linear Fit
Xinv = pinv([ones(prod(size(x)),1),x(:)]);
betaLinear = Xinv * yD(:);

figure

plot(x,yD','o','Color',[0.5 0.5 0.5])
hold all
xPlot = linspace(0.5,10.5,20);
plot(xPlot,xPlot*betaLinear(2)+betaLinear(1),'k','LineWidth',3)
xlabel('predictor value (e.g. saccadic amplitude)')
ylabel('signal to be explained')

set(get(gcf,'Children'),'Box','off','YTick',[],'XTick',[])

set(gcf,'Position',[   716   687   763   309])
export_fig nonlinear1.png -transparent
%% Boxcar function fit
tmp = floor(x);
basisBox = [];
%build the dummy coding matrix
for u = unique(tmp)';
basisBox = [basisBox tmp == u];
end
betaBox = pinv(basisBox*1)*yD(:);


figure,
subplot(1,2,1)
plot(x,basisBox','-')
ylim([0 3])
xlim([0 11])
xlabel('basis set')

subplot(1,2,2)
plot(x,yD','o','Color',[0.5 0.5 0.5])

hold all

plot(linspace(min(x),max(x),N),basisBox*betaBox,'k','LineWidth',3)
set(get(gcf,'Children'),'Box','off','YTick',[],'XTick',[])

xlabel('predictor value (e.g. saccadic amplitude)')
ylabel('signal to be explained')
set(gcf,'Position',[   716   687   763   309])
export_fig nonlinear_boxcar.png -transparent
%% Spline Fit
tmp = ceil(x*19/7)-2;

knotseq = [min(x) min(x) min(x) linspace(min(x),max(x),4) max(x) max(x) max(x)];
splinePlot = Bernstein(linspace(min(x),max(x),200),knotseq,[],4);
splineFit= Bernstein(x,knotseq,[],4);
splineFit(1,1) = 1;
splineFit(end,end)= 1;
betaSpline = pinv([splineFit]) * yD(:)

figure
subplot(1,2,1)
plot(splinePlot,'-')
xlabel('basis set')
subplot(1,2,2)
plot(x,yD','o','Color',[0.5 0.5 0.5])
hold all,
plot(0,0),plot(0,0),plot(0,0),plot(0,0),plot(0,0),plot(0,0); %stupid colors were not the same in both plots
plot(linspace(min(x),max(x),200),splinePlot'.*betaSpline(1:end),'-')


plot(linspace(min(x),max(x),200),splinePlot*betaSpline(1:end),'k','LineWidth',3)

xlabel('predictor value (e.g. saccadic amplitude)')
ylabel('signal to be explained')

set(get(gcf,'Children'),'Box','off','YTick',[],'XTick',[])
set(gcf,'Position',[   716   687   763   309])
export_fig nonlinear_spline.png -transparent