clear all
clc
 
sigma = 0.01;

%truncationError = zeros(1,optionalNum);
xNum  = 200;
% xNum  = 301;
tNum  = floor(xNum^(7/8));
xMax  = 1; % xMin = 0
tMax  = 2;
dx    = xMax/(xNum-1);
dt    = tMax/(tNum-1);
xData = 0:dx:xMax;
tData = 0:dt:tMax;
[xMesh,tMesh] = meshgrid(xData,tData); % col dim: x; row dim: t

% True solution: F(x-ct)
 
c = 2.5; % Specify the convection speed
F  = @(x) 2*sin(pi*x/4); % Specify the function F
dF = @(x) cos(pi*x/4)*pi/2;
u   = F(xMesh-c*tMesh);
rng(3)
uNoise = u + normrnd(0,sigma,size(u));

% tic
%  uDenoise = LocalPolyRegression(uNoise,2,1*tNum^(-1/7),1.2*xNum^(-1/8),tMesh,xMesh,3);
% toc


tic
uDenoise = FastLocalPolyRegression(uNoise,2,1*tNum^(-1/7),0.8*xNum^(-1/8),dt,dx,3);
toc

diffdUx = (uNoise(:,2:end)-uNoise(:,1:end-1))/dx;
diffdUt = (uNoise(2:end,:)-uNoise(1:end-1,:))/dt;
% uDenoise = LocalPolyRegression(uNoise,2,tNum^(-1/7),xNum^(-1/8),tMesh,xMesh,3);
B   = 0;
denoisedUt = uDenoise(:,B+1:end-B,end); 
denoiseU = uDenoise(:,B+1:end-B,1); 
denoiseU2 = denoiseU.*denoiseU;
denoisedUx = uDenoise(:,B+1:end-B,2);
denoisedUx2 = denoisedUx.*denoisedUx;
denoiseUdUx = denoisedUx.*denoiseU;
denoisedUxx = uDenoise(:,B+1:end-B,3);
denoiseUdUxx = denoisedUxx.*denoiseU;
denoisedUxdUxx = denoisedUxx.*denoisedUx;
denoisedUxx2 = denoisedUxx.*denoisedUxx;
constTerm = ones(size(denoiseU));

dUtVec   = denoisedUt(:);
constVec = constTerm(:);
UVec = denoiseU(:);
U2Vec = denoiseU2(:);
dUxVec = denoisedUx(:);
dUx2Vec = denoisedUx2(:);
UdUxVec = denoiseUdUx(:);
dUxxVec = denoisedUxx(:);
dUxx2Vec = denoisedUxx2(:);
dUxdUxxVec = denoisedUxdUxx(:);
UdUxxVec   = denoiseUdUxx(:);

%% u vs uhat
figure
tn = ceil(1/dt);
plot(xData,u(tn,:),'k-.','LineWidth',1.5)
hold on
plot(xData,uNoise(tn,:),'r.','MarkerSize',10,'Color',[0.8500 0.3250 0.0980])
plot(xData,denoiseU(tn,:),'-','LineWidth',1.5,'Color',[0 0.4470 0.7410])
hleg1 = legend({'$u(\cdot,1)$', 'noisy data',  '$\widehat{u}(\cdot,1)$'},'Interpreter','Latex','Location','northwest');
xlabel('$x$','Interpreter','Latex')
title({'$u$ Estimation ($\sigma=0.01$)'},'Interpreter','Latex')
set(gca,'TickLabelInterpreter', 'Latex');
set(gca,'fontsize',30)
%% dux vs duhatx
dux = dF(xMesh-c*tMesh);
figure
yyaxis left
h1=plot(xData,dux(tn,:),'k-.','LineWidth',1.5);
hold on
yyaxis right
h2=plot(xData(1:end-1),diffdUx(tn,:),'.','MarkerSize',10);
ylim([-10,10])
yyaxis left
h3=plot(xData,denoisedUx(tn,:),'-','LineWidth',1.5,'Color',[0 0.4470 0.7410]);
hleg1 = legend([h1,h2,h3],{'$\partial_xu(\cdot,1)$', 'Forward Diff',  '$\widehat{\partial_xu}(\cdot,1)$'},'Interpreter','Latex','Location','northwest');
xlabel('$x$','Interpreter','Latex')
title({'$\partial_xu$ Estimation ($\sigma=0.01$)'},'Interpreter','Latex')
set(gca,'TickLabelInterpreter', 'Latex');
set(gca,'fontsize',30)
%% dut vs duhatt
dut = dF(xMesh-c*tMesh)*(-c);
figure
xn = ceil(0.5/dx);
yyaxis left
h1=plot(tData,dut(:,xn),'k-.','LineWidth',1.5);
ylim([-6,6])
hold on
yyaxis right
h2=plot(tData(1:end-1),diffdUt(:,xn),'.','MarkerSize',10);
ylim([-6,6])
yyaxis left
h3=plot(tData,denoisedUt(:,xn),'-','LineWidth',1.5,'Color',[0 0.4470 0.7410]);
hleg1 = legend([h1,h2,h3],{'$\partial_tu(0.5,\cdot)$', 'Forward Diff',  '$\widehat{\partial_tu}(0.5,\cdot)$'},'Interpreter','Latex','Location','northwest');
xlabel('$t$','Interpreter','Latex')
title({'$\partial_tu$ Estimation ($\sigma=0.01$)'},'Interpreter','Latex')
set(gca,'TickLabelInterpreter', 'Latex');
set(gca,'fontsize',30)
%%
Nspace = linspace(10,10^8,10^7);
figure
sigma = 0.3;
for K = [10,50,100,200]
    PP = (8*K+2).*Nspace.^(15/7).*K.*exp(-(Nspace.^(1/7)-3).^2/(2*sigma^2));
    loglog(Nspace,PP,'LineWidth',1.5)
    hold on
end
line([10,10^8],[1,1],'LineStyle','-.','Color',[0.5,0.5,0.5],'LineWidth',1.5)
xlim([10,10^8])
ylim([10^-20,10^20])
xlabel('N','Interpreter','Latex')

set(gca,'fontsize',20)
%%
clear all
clc
 
sigma = 0.5;

%truncationError = zeros(1,optionalNum);
xNumList = [300];
optionalNum = length(xNumList);
NormList = zeros(1,optionalNum);
k=1;
xNum  = ceil(xNumList(k));
% xNum  = 301;
tNum  = floor(xNum^(7/8));
xMax  = 10; % xMin = 0
tMax  = 10;
dx    = xMax/(xNum-1);
dt    = tMax/(tNum-1);
xData = 3*dx:dx:xMax;
tData = dt:dt:tMax;
[xMesh,tMesh] = meshgrid(xData,tData); % col dim: x; row dim: t
Xmax = 0.5;
Tmax = 0.4;

A = 0.1;
ratio_constant = (8./(xMesh.*A)).*abs(cos(xMesh*Xmax-tMesh*Tmax)-cos(tMesh*Tmax)-cos(xMesh*Xmax)+1)./...
    ((4*xMesh.*tMesh*Xmax*Tmax+1+cos(2*xMesh*Xmax-2*tMesh*Tmax)-cos(2*xMesh*Xmax)-cos(2*tMesh*Tmax)));

ratio_u = (1./(xMesh)).*abs(sin(2*xMesh*Xmax-2*tMesh*Tmax)+sin(2*tMesh*Tmax)-sin(2*xMesh*Xmax))./...
    ((4*xMesh.*tMesh*Xmax*Tmax+1+cos(2*xMesh*Xmax-2*tMesh*Tmax)-cos(2*xMesh*Xmax)-cos(2*tMesh*Tmax)));

ratio_uxx = xMesh.*abs(sin(2*xMesh*Xmax-2*tMesh*Tmax)+sin(2*tMesh*Tmax)-sin(2*xMesh*Xmax))./...
    ((4*xMesh.*tMesh*Xmax*Tmax+1+cos(2*xMesh*Xmax-2*tMesh*Tmax)-cos(2*xMesh*Xmax)-cos(2*tMesh*Tmax)));
%%
figure
[C,h] = contourf(tMesh,xMesh,max(max(ratio_u,ratio_constant),ratio_uxx),[0,0.5,1],'LineWidth',1.5,'ShowText','on');
colormap(summer)
clabel(C,h,'FontSize',20,'Color','k','FontWeight','bold','Interpreter','Latex')
% legend({'$H(\omega,c) <1$'},'Interpreter','Latex');
% ylim([3*dx,0.5])
hold on 
xlabel('$c$','interpreter','Latex')
ylabel('$\omega$','interpreter','Latex')
set(gca,'fontsize',20)