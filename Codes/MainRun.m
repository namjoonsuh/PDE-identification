clear all
clc
sigma = 0.01;
xNum  = 100;
tNum  = floor(xNum^(7/8));
xMax  = 1; 
tMax  = 0.1;
dx    = xMax/xNum;
dt    = tMax/tNum;
xData = 0:dx:xMax;
tData = 0:dt:tMax;
[xMesh,tMesh] = meshgrid(xData,tData); % col dim: x; row dim: t

fineRatioX = 1; 
fineRatioT = 5000;
fineDx = dx/fineRatioX;
fineDt = dt/fineRatioT;
fineXData = 0:fineDx:xMax;
fineTData = 0:fineDt:tMax;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%% Run Parabolic Example: Heat equation (0-Dirichlet boundary)
% % f: initial condition
% % nu: conductivity parameter
% % FTCS scheme
% f = @(x) 5*cos(5*pi*x).*(1-x).^2.*x;
% nu = 0.1;
% U = ParabolicGen(f,fineDx,fineDt,xMax,tMax,nu);

%%% Run Hyperbolic Example
% f: initial condition
% a: convection parameter
% Lax-Wendroff scheme
f = @(x) 3*sin(4*pi*x);
a = 1;
U = HyperbolicGen(f,fineDx,fineDt,xMax,tMax,a);

% %%% Run Advection-Diffusion Example
% f = @(x) sin(4*pi*x);
% a = 2;
% U = ADDGen(f,fineDx,fineDt,xMax,tMax,a,0.1);

% %%% Run Bateman-Burgers Example
% f = @(x) sin(4*pi*x).^2+sin(2*pi*x).^3;
% % f = @(x) sin(2*pi*x);
% nu = 0.01;
% U = BurgersGen(f,fineDx,fineDt,xMax,tMax,nu);

% %%% Linearized KdV Example
% f = @(x) cos(pi*x);
% U = KdVGen(f,fineDx,fineDt,xMax,tMax);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Downsample
if fineRatioX==1
    U = U(mod(1:length(fineTData)-1,fineRatioT)==1,:);
else
    U = U(mod(1:length(fineTData)-1,fineRatioT)==1,...
        mod(1:length(fineXData),fineRatioX)==1);
end 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(3)
uNoise = U + normrnd(0,sigma,size(U));
uDenoise = FastLocalPolyRegression(uNoise,2,0.15*tNum^(-1/7),0.5*xNum^(-1/8),dt,dx,3);

diffdUx = (uNoise(:,2:end)-uNoise(:,1:end-1))/dx;
diffdUt = (uNoise(2:end,:)-uNoise(1:end-1,:))/dt;
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

evalMat =  max(abs((noiseFmat(:,[1:3,5:end])'*dUxVec)*(dUxVec'*dUxVec)^(-1)));
[LL,fitinfo] = lasso(noiseFmat,dUtVec);

%%

figure
for i = 1:size(LL,1)
    if i == 6
        h1 = semilogx(fitinfo.Lambda,LL(i,:),'LineWidth',2,'Color',[0.8500 0.3250 0.0980]);
    elseif i == 7
        h2 = semilogx(fitinfo.Lambda,LL(i,:),'LineWidth',2,'Color',[0.4660 0.6740 0.1880]);
    else    
        semilogx(fitinfo.Lambda,LL(i,:),'-.','Color',[0.5,0.5,0.5],'LineWidth',1.5)
    end
    hold on
end
xlim([0,2])
% ylim([-1.1,0.1])

xlabel('$\lambda$','Interpreter','Latex')
ylabel('Coefficient','Interpreter','Latex')
hleg1 = legend([h1],{'$uu_{x}$'},'Interpreter','Latex');
set(gca,'TickLabelInterpreter','latex')

% title('$\sigma=0.03,M=1000$','Interpreter','Latex')

set(gca,'fontsize',30)
%%
figure
for i = floor(linspace(1,tNum,5))
    if i == 1
        plot(xData,uNoise(i,:),'Color',[0,0,0],'LineWidth',1.5)
    else
        plot(xData,uNoise(i,:),'Color',[0 0.4470 0.7410,i/tNum],'LineWidth',1.5)
    end
    hold on 
end
ylim([-2,2.5])
xlabel('$x$','Interpreter','Latex')
ylabel('$\widetilde{u}(x,\cdot)$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',30)