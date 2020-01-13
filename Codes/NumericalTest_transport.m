clear all
clc
 
sigma = 0.01;

%truncationError = zeros(1,optionalNum);
xNumList = [10,50,100,200,300,400,500,600,700,800,900,1000];
optionalNum = length(xNumList);
NormList = zeros(100,optionalNum);
SuccessList = zeros(100,optionalNum);
for k = 1:optionalNum
for experiment = 1:100
xNum  = ceil(xNumList(k));
% xNum  = 301;
tNum  = floor(xNum^(7/8));
xMax  = 1; % xMin = 0
tMax  = 1;
dx    = xMax/(xNum-1);
dt    = tMax/(tNum-1);
xData = 0:dx:xMax;
tData = 0:dt:tMax;
[xMesh,tMesh] = meshgrid(xData,tData); % col dim: x; row dim: t

% True solution: F(x-ct)
 
c = 2.5/2; % Specify the convection speed
% F  = @(x) 2*sin(pi*x/4); % Specify the function F
% dF = @(x) cos(pi*x/4)*pi/2;
% u   = F(xMesh-c*tMesh);
u = 2*sin(4*xMesh-8*tMesh);
% rng(3)
uNoise = u + normrnd(0,sigma,size(u));
% tic
% uDenoise = LocalPolyRegression(uNoise,2,2*tNum^(-1/7),3*xNum^(-1/8),tMesh,xMesh,3);
% toc 
% tic
uDenoise = FastLocalPolyRegression(uNoise,2,1*tNum^(-1/7),1.5*xNum^(-1/8),dt,dx,3);
% toc
%%
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

noiseFmat = [constVec,UVec,U2Vec,...
    dUxVec,dUx2Vec,UdUxVec,dUxxVec,...
    dUxx2Vec,dUxdUxxVec,UdUxxVec];

% noiseFmat = [constVec,UVec,...
%     dUxVec,dUxxVec];

evalMat =  (noiseFmat(:,[1:3,5:end])'*dUxVec)*(dUxVec'*dUxVec)^(-1);
NormList(experiment,k) = max(sum(abs(evalMat),2));
[LL,fitinfo] = lasso(noiseFmat,dUtVec);
SuccessList(experiment,k) = (find(LL(:,end-1)~=0)==4);
%truncationError(k) = max(abs(dUtVec+2.5*dUxVec));
% figure
% surf(xData,tData,reshape(abs(dUtVec+2.5*dUxVec),tNum,xNum))
% zlim([0,0.5])
disp(sprintf('experiment %d for xNum = %d', experiment,xNumList(k)))
end
end
%% MIP test plot
figure
xNumList = [10,50,100,200,300,400,500,600,700,800,900,1000];
load('transportMIP_0_3.mat')
plot1 = plot(xNumList,NormList,'LineWidth',1,'Color',[0.5,0.5,0.5,0.1]);
xlim([10,1000])
hold on
h = boxplot(NormList, 'positions', xNumList,'labels', xNumList);
set(h,{'linew'},{2})
set(gca,'YScale','log')
xlabel('M','Interpreter','Latex')
plot2 = line([10,1000],[1,1],'LineStyle','-.','LineWidth',1,'Color','r');
alpha(plot2,1)
load('transportMIP_0.mat')
plot3 = plot(xNumList,NormList,'LineStyle','-.','LineWidth',1.5,'Color',[0.4940 0.1840 0.5560]);
ylabel('$||\widehat{\mathbf{F}}_{\mathcal{S}^c}^T\widehat{\mathbf{F}}_{\mathcal{S}}\big(\widehat{\mathbf{F}}_{\mathcal{S}}^T\widehat{\mathbf{F}}_{\mathcal{S}}\big)^{-1}||_\infty$','Interpreter','Latex')
% line([0,4000],[1,1],'Color',[0.5,0.5,0.5],'LineWidth',1,'LineStyle','-.')
hleg1 = legend([plot2,plot3],{'$\ell_\infty$-norm $=1$','$\sigma=0$'},'Interpreter','Latex');
% xlim([20,400])
% title('MIP for Transport ($\sigma=0.3$)','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',20)
set(gcf, 'Position',  [100, 100, 800, 300])
ylim([0,10^2])
%% u vs uhat
figure
tn = 30;
plot(xData,u(tn,:),'k-.','LineWidth',1.5)
hold on
plot(xData,uNoise(tn,:),'.','MarkerSize',10)
plot(xData,denoiseU(tn,:),'-','LineWidth',1.5)
hleg1 = legend({'$u(\cdot,30\Delta t)$', 'noisy data',  '$\widehat{u}(\cdot,30\Delta t)$'},'Interpreter','Latex','Location','northwest');
xlabel('$x$','Interpreter','Latex')
title({'Estimating $u$ using Noisy Data ($\sigma=0.01$)'},'Interpreter','Latex')
set(gca,'fontsize',20)
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
h3=plot(xData,denoisedUx(tn,:),'-','LineWidth',1.5);
hleg1 = legend([h1,h2,h3],{'$\partial_xu(\cdot,30\Delta t)$', 'Forward Diff',  '$\widehat{\partial_xu}(\cdot,30\Delta t)$'},'Interpreter','Latex','Location','northwest');
xlabel('$x$','Interpreter','Latex')
title({'Estimating $\partial_xu$ using Noisy Data ($\sigma=0.01$)'},'Interpreter','Latex')
set(gca,'fontsize',20)
%% dut vs duhatt
dut = dF(xMesh-c*tMesh)*(-c);
figure
yyaxis left
h1=plot(tData,dut(:,tn),'k-.','LineWidth',1.5);
ylim([-6,6])
hold on
yyaxis right
h2=plot(tData(1:end-1),diffdUt(:,tn),'.','MarkerSize',10);
ylim([-6,6])
yyaxis left
h3=plot(tData,denoisedUt(:,tn),'-','LineWidth',1.5);
hleg1 = legend([h1,h2,h3],{'$\partial_tu(30\Delta x,\cdot)$', 'Forward Diff',  '$\widehat{\partial_tu}(30\Delta x,\cdot)$'},'Interpreter','Latex','Location','northwest');
xlabel('$t$','Interpreter','Latex')
title({'Estimating $\partial_tu$ using Noisy Data ($\sigma=0.01$)'},'Interpreter','Latex')
set(gca,'fontsize',20)
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