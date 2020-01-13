clear all
clc
 
sigma = 1;

%truncationError = zeros(1,optionalNum);
xNumList = [50,100,200,300,400,500,600,700,800,900,1000,1200,1500,2000,2500,3000];
optionalNum = length(xNumList);
lambdaList = zeros(1,optionalNum);
for k = 1:optionalNum
xNum  = ceil(xNumList(k));
% xNum  = 301;
% tNum  = 31;
tNum  = floor(xNum^(7/8));
xMax  = 2; % xMin = 0
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
u = 2*sin(3*xMesh-8*tMesh);
rng(3)
uNoise = u + normrnd(0,sigma,size(u));
uDenoise = FastLocalPolyRegression(uNoise,2,0.15*tNum^(-1/7),3*xNum^(-1/8),dt,dx,3);
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


[LL,fitinfo] = lasso(noiseFmat,dUtVec,'Lambda',linspace(0,2.5,100));

llist  = find((sum(LL~=0)==1)==1);
lambdaList(k) = fitinfo.Lambda(llist(1))
find(LL(:,llist(1))~=0)

end
%%
load('lambdaListAdv_1.mat')
figure
tNumList = xNumList.^(7/8);
h1 = plot(tNumList(1:end),lambdaList(1:end),'k','LineWidth',1.5);
hold on 
ratio =  lambdaList(1:end)*(sqrt(log(tNumList(1:end))./tNumList(1:end).^(4/7)))'/norm(sqrt(log(tNumList(1:end))./tNumList(1:end).^(4/7)))^2;
h2 = plot(tNumList(1:end),sqrt(log(tNumList(1:end))./tNumList(1:end).^(4/7))*ratio,'r','LineWidth',1.5);
xlim([tNumList(1),tNumList(end)])
hleg1 = legend([h1,h2],{'$\lambda_{\min}$','$\sim\sqrt{\frac{\ln N}{N^{4/7}}}$'},'Interpreter','Latex');
xlabel('N','Interpreter','Latex')
set(gca,'fontsize',20,'TickLabelInterpreter','Latex')





