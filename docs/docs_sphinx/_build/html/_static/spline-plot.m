b = [0.1 0.5 1 2 0.2 4 2]/2;
N = 200;
x = sort(rand(N,1)*7+1);
y = [];
yD = log(8-x) + randn(N,1)*0.1;
% yD = yD-min(yD);

%% Linear Fit
Xinv = pinv([ones(prod(size(x)),1),x(:)]);
betaLinear = Xinv * yD(:);

figure

plot(x,yD','o','Color',[0.5 0.5 0.5])
hold all
xPlot = linspace(1,11,20);
plot(xPlot,xPlot*betaLinear(2)+betaLinear(1),'k','LineWidth',3)
xlabel('predictor value (e.g. saccadic amplitude)')
ylabel('signal to be explained')

set(get(gcf,'Children'),'Box','off','YTick',[],'XTick',[])

set(gcf,'Position',[   716   687   763   309])
export_fig nonlinear1.png -transparent
%% Boxcar function fit
tmp = floor(x)

%build the dummy coding matrix
basisBox = [tmp == 1 tmp==2 tmp==3 tmp==4 tmp==5 tmp==6 tmp==7];
betaBox = pinv(basisBox*1)*yD(:)


figure,
subplot(1,2,1)
plot(x,basisBox','-')
ylim([0 3])
xlim([1 8])
xlabel('basis set')

subplot(1,2,2)
plot(x,yD','o','Color',[0.5 0.5 0.5])

hold all
xPlot = linspace(1.1,8,100)';
yTmp = [];
for k = 2:1:8.5
   yTmp =[yTmp   ((k-1) < xPlot) & (xPlot <= k)];
end

plot(xPlot,yTmp*betaBox,'k','LineWidth',3)
set(get(gcf,'Children'),'Box','off','YTick',[],'XTick',[])

xlabel('predictor value (e.g. saccadic amplitude)')
ylabel('signal to be explained')
set(gcf,'Position',[   716   687   763   309])
export_fig nonlinear_boxcar.png -transparent
%% Spline Fit
basisSpline = Bernstein(1:20,[1 1 1 linspace(1,20,5) 20 20 20],[],4);
tmp = ceil(x*19/7)-2;
basisBox = [tmp == 1 tmp==2 tmp==3 tmp==4 tmp==5 tmp==6 tmp==7 tmp==8 tmp==9 tmp==10 tmp==11 tmp==12 tmp==13 tmp==14 tmp==15 tmp==16 tmp == 17 tmp==18 tmp==19 tmp ==20];

betaSpline = pinv([ones(N,1) basisBox*basisSpline]) * yD(:)

figure
subplot(1,2,1)
plot(basisSpline,'--')
xlabel('basis set')
subplot(1,2,2)
plot(x,yD','o','Color',[0.5 0.5 0.5])
hold all,
plot(linspace(1,8,20),basisSpline'.*betaSpline(2:end),'--')
plot(linspace(1,8,20),basisSpline*betaSpline(2:end)+betaSpline(1),'k','LineWidth',3)

xlabel('predictor value (e.g. saccadic amplitude)')
ylabel('signal to be explained')

set(get(gcf,'Children'),'Box','off','YTick',[],'XTick',[])
set(gcf,'Position',[   716   687   763   309])
export_fig nonlinear_spline.png -transparent