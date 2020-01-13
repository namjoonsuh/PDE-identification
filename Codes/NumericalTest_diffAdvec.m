clear all
clc
 
% tic
sigma = 1;

%truncationError = zeros(1,optionalNum);
xNumList = [50,100,200,300,400,500,600,700,800,900,1000];
optionalNum = length(xNumList);
NormList = zeros(100,optionalNum);
for k = 1:optionalNum
for experiment = 1:100

xNum  = ceil(xNumList(k));
% xNum  = 301;
tNum  = floor(xNum^(7/8));
xMax  = 3; % xMin = 0
tMax  = 0.5;
dx    = xMax/(xNum-1);
dt    = tMax/(tNum-1);
xData = 0:dx:xMax;
tData = 0:dt:tMax;
[xMesh,tMesh] = meshgrid(xData,tData); % col dim: x; row dim: t

% True solution: sqrt(a)/sqrt(a+4*D*t)*exp(-(x-v*t)^2/(a+4*D*t))
a = 0.6;
D = 2; % Specify the diffusitivity
v = 10; % Specify the convection speed
u = sqrt(a)./sqrt(a+4*D*tMesh).*exp(-((xMesh+0.7)-v*tMesh).^2./(a+4*D*tMesh));

uNoise = u + normrnd(0,sigma,size(u));
uDenoise = FastLocalPolyRegression(uNoise,2,0.2*tNum^(-1/7),4*xNum^(-1/8),dt,dx,3);
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

evalMat =  (noiseFmat(:,[1:3,5,6,8:end])'*noiseFmat(:,[4,7]))*(noiseFmat(:,[4,7])'*noiseFmat(:,[4,7]))^(-1);
NormList(experiment,k) = max(sum(abs(evalMat),2));
disp(sprintf('experiment %d for xNum = %d', experiment,xNumList(k)))
%truncationError(k) = max(abs(dUtVec+2.5*dUxVec));
% figure
% surf(xData,tData,reshape(abs(dUtVec+2.5*dUxVec),tNum,xNum))
% zlim([0,0.5])
end
end


%% MIP test plot
figure
xNumList = [50,100,200,300,400,500,600,700,800,900,1000,1500,2000];
load('advecdiff_0_1.mat')
plot1 = plot(xNumList,NormList,'LineWidth',1,'Color',[0.5,0.5,0.5,0.1]);
xlim([10,1000])
hold on
h = boxplot(NormList, 'positions', xNumList,'labels', xNumList);
set(h,{'linew'},{2})
set(gca,'YScale','log')
xlabel('M','Interpreter','Latex')
plot2 = line([10,1000],[1,1],'LineStyle','-.','LineWidth',1,'Color','r');
alpha(plot2,1)
load('advecdiff_0.mat')
plot3 = plot(xNumList,NormList,'LineStyle','-.','LineWidth',1.5,'Color',[0.4940 0.1840 0.5560]);
ylabel('$||\widehat{\mathbf{F}}_{\mathcal{S}^c}^T\widehat{\mathbf{F}}_{\mathcal{S}}\big(\widehat{\mathbf{F}}_{\mathcal{S}}^T\widehat{\mathbf{F}}_{\mathcal{S}}\big)^{-1}||_\infty$','Interpreter','Latex')
% line([0,4000],[1,1],'Color',[0.5,0.5,0.5],'LineWidth',1,'LineStyle','-.')
hleg1 = legend([plot2,plot3],{'$\ell_\infty$-norm $=1$','$\sigma=0$, w. LPR, $\sigma=0$ wo. LPR'},'Interpreter','Latex');
% xlim([20,400])
title('MIP for Transport ($\sigma=0.3$)','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',20)
set(gcf, 'Position',  [100, 100, 800, 300])
ylim([0,10^2])


% toc
% figure
% plot1 = plot(xNumList,NormList,'LineWidth',1.5);