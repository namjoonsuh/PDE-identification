%%%%%%%%%% Generate Hyperbolic Example Data %%%%%%%%%%%%%%%%
%%%%%%%%%% Roy  Jan.2.2020 %%%%%%%%
%%%% Periodic boundary condition
%%%% Lax-Wendroff Scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = HyperbolicGen(f,dx,dt,xMax,tMax,a)
M = floor(xMax/dx);
N = floor(tMax/dt);
U = zeros(N+1,M+1);
U(1,:) = f(0:dx:dx*M);
R = a*dt/dx;
if R > 1
    warning(sprintf('Your choice of grid leads to unstable scheme: R=%4f',R))
end
for n = 2:N+1
    U(n,2:end-1) = U(n-1,2:end-1)-...
        R/2*(U(n-1,3:end)-U(n-1,1:end-2))+...
        R^2/2*(U(n-1,3:end)-2*U(n-1,2:end-1)+U(n-1,1:end-2));
    U(n,1) = U(n-1,1)-...
        R/2*(U(n-1,2)-U(n-1,end-1))+...
        R^2/2*(U(n-1,2)-2*U(n-1,1)+U(n-1,end-1));
    U(n,end) = U(n,1);
end

end