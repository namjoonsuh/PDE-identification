%%%%%%%%%% Generate Fisher's Example Data %%%%%%%%%%%%%%%%
%%%%%%%%%% Roy  Jan.2.2020 %%%%%%%%
%%%% Periodic boundary condition
%%%% Zabusky Kruskal scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = FisherGen(f,dx,dt,xMax,tMax,D,r)
M = floor(xMax/dx);
N = floor(tMax/dt);

U = zeros(N+1,M+1);
U(1,:) = f(0:dx:dx*M);
R = D*dt/dx^2;
for n = 2:N+1
    U(n,2:end-1) = U(n-1,2:end-1)+...
        R*(U(n-1,3:end)-2*U(n-1,2:end-1)+U(n-1,1:end-2))+...
        r*U(n-1,2:end-1).*(1-U(n-1,2:end-1))*dt;
end
end