%%%%%%%%%% Local Polynomial Regression %%%%%%%%%%%%%%%%
%%%%%%%%%% Roy  Dec.1.2019 %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function smoothedData = LocalPolyRegression(u,maxOrder,h,w,tMesh,xMesh,D)
%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%
% u       : function value timeDim * spaceDim
% maxOrder: maximum order of the spatial derivatives 
% h       : bandwidth for time  
% w       : bandwidth for space
% tMesh   : time mesh
% xMesh   : space mesh
% D       : extra order needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K     = @(x) subplus(3/4*(1-((xMesh(1,:)-x)/w).^2));
%K     = @(x) subplus(1-abs((xMesh(1,:)-x)/w));
smoothedData = zeros(size(u,1),size(u,2),maxOrder+2);
for order = 1:maxOrder+1
    for n = 1:size(u,1)
        Y = u(n,:)';
        for i = 1:size(u,2)
            W = K(xMesh(n,i))';
            X = ((xMesh(1,:)-xMesh(1,i))').^(0:(order+D));
            b = (X'*(W.*X))\(X'*(W.*Y));
            smoothedData(n,i,order) = factorial(order-1)*b(order);
        end
    end
end

K     = @(t) subplus(3/4*(1-((tMesh(:,1)-t)/h).^2));
%K     = @(t) subplus(1-abs((tMesh(:,1)-t)/h));
for m = 1:size(u,2)
    Y = u(:,m);
    for i = 1:size(u,1)
        W = K(tMesh(i,m));
        X = (tMesh(:,1)-tMesh(i,1)).^(0:4);
        b = (X'*(W.*X))\(X'*(W.*Y));
        smoothedData(i,m,maxOrder+2) = b(2);
    end
end

end