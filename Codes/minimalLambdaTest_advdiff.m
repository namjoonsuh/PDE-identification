clear all
clc
 
sigma = 0.05;

%truncationError = zeros(1,optionalNum);
xNumList = [50,100,200,300,400,500,600,700,800,900,1000,1200,1500,2000,2500,3000];
optionalNum = length(xNumList);
lambdaList = zeros(1,optionalNum);
for k = 1:optionalNum
xNum  = ceil(xNumList(k));
% xNum  = 301;
tNum  = floor(xNum^(7/8));
xMax  = 3; % xMin = 0
tMax  = 0.4;
dx    = xMax/(xNum-1);
dt    = tMax/(tNum-1);
xData = 0:dx:xMax;
tData = 0:dt:tMax;
[xMesh,tMesh] = meshgrid(xData,tData); % col dim: x; row dim: t

% True solution: sqrt(a)/sqrt(a+4*D*t)*exp(-(x-v*t)^2/(a+4*D*t))
a = 0.9;
D = 5; % Specify the diffusitivity
v = 11; % Specify the convection speed
u = sqrt(a)./sqrt(a+4*D*tMesh).*exp(-((xMesh+1.6)-v*tMesh).^2./(a+4*D*tMesh));
rng(1)
uNoise = u + normrnd(0,sigma,size(u));
uDenoise = FastLocalPolyRegression(uNoise,2,0.15*tNum^(-1/7),3*xNum^(-1/8),dt,dx,3);
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



[LL,fitinfo] = lasso(noiseFmat,dUtVec,'Lambda',linspace(0,2.5,100));
llist = find(sum((LL~=0).*[1 1 1 0 1 1 0 1 1 1]')==0);
lambdaList(k) = fitinfo.Lambda(llist(1));
find(LL(:,llist(1))~=0)

% evalMat =  (noiseFmat(:,[1:3,5,6,8:end])'*noiseFmat(:,[4,7]))*(noiseFmat(:,[4,7])'*noiseFmat(:,[4,7]))^(-1);
% NormList(k) = max(sum(abs(evalMat),2));
% %truncationError(k) = max(abs(dUtVec+2.5*dUxVec));
% % figure
% % surf(xData,tData,reshape(abs(dUtVec+2.5*dUxVec),tNum,xNum))
% % zlim([0,0.5])
end

%%
tNumList = xNumList.^(7/8);
h1 = plot(tNumList(3:end),lambdaList(3:end),'k','LineWidth',1.5);
hold on 
h2 = plot(tNumList(3:end),sqrt(log(tNumList(3:end))./tNumList(3:end).^(4/7))/1.15,'r','LineWidth',1.5);
xlim([tNumList(3),tNumList(end)])
hleg1 = legend([h1,h2],{'$\lambda_{\min}$','$\sim\sqrt{\frac{\ln N}{N^{4/7}}}$'},'Interpreter','Latex');
xlabel('N','Interpreter','Latex')
set(gca,'fontsize',20,'TickLabelInterpreter','Latex')