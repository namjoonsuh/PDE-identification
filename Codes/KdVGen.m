%%%%%%%%%% Generate KdV Example Data %%%%%%%%%%%%%%%%
%%%%%%%%%% Roy  Jan.2.2020 %%%%%%%%
%%%% Periodic boundary condition
%%%% Zabusky Kruskal scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = KdVGen(f,dx,dt,xMax,tMax)
M = floor(xMax/dx);
N = floor(tMax/dt);

U = zeros(N+1,M+1);
U(1,:) = f(0:dx:dx*M);
R = dt/dx;
r = dt/dx^3;
for n = 2:N+1
    U(n,:) = U(n-1,:)-r/2*(U(n-1,[3:end,1:2])-2*U(n-1,[2:end,1])+...
        2*U(n-1,[end,1:end-1])-U(n-1,[end-1,end,1:end-2]));
end
end