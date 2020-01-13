%%%%%%%%%% Fast Local Polynomial Regression %%%%%%%%%%%%%%%%
%%%%%%%%%% Roy  Jan.2.2020 %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function smoothedData = FastLocalPolyRegression(u,maxOrder,h,w,dt,dx,D)
%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%
% u       : function value timeDim * spaceDim
% maxOrder: maximum order of the spatial derivatives 
% h       : bandwidth for time  
% w       : bandwidth for space
% tMesh   : time mesh
% xMesh   : space mesh
% D       : extra order needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K     = @(x) subplus(3/4*(1-(x/w).^2));
%K     = @(x) subplus(1-abs((xMesh(1,:)-x)/w));
smoothedData = zeros(size(u,1),size(u,2),maxOrder+2);
halfWindowNum = ceil(w/dx);
gridList = 1:halfWindowNum;
gridList = [-fliplr(gridList),0,gridList]*dx;
weight = K(gridList);
X = fliplr(vander(gridList));
X = X(:,1:maxOrder+D);

for n = 1:size(u,1)
    Yextend  = padarray(u(n,:),[0,halfWindowNum],'symmetric');
    Y = circshift(flipud(gallery('circul',Yextend)),1);
    Y = Y(1:2*halfWindowNum+1,1:size(u,2));
    beta = X'*(weight'.*Y);
    beta = (X'*(weight'.*X))\beta;
    smoothedData(n,:,1:end-1) = (factorial(0:maxOrder)'.*beta(1:(maxOrder+1),:))';
end



K     = @(x) subplus(3/4*(1-(x/h).^2));
halfWindowNum = ceil(h/dt);
gridList = 1:halfWindowNum;
gridList = [-fliplr(gridList),0,gridList]*dt;

weight = K(gridList);
X = fliplr(vander(gridList));
if 1+D<size(X,2)
    X = X(:,1:1+D);
end


for m = 1:size(u,2)
    Yextend  = padarray(u(:,m),[halfWindowNum,0],'symmetric');
    Y = circshift(flipud(gallery('circul',Yextend)),1);
    Y = Y(1:2*halfWindowNum+1,1:size(u,1));
    beta = X'*(weight'.*Y);
    beta = (X'*(weight'.*X))\beta;
    smoothedData(:,m,end) = beta(2,:)';
end

end