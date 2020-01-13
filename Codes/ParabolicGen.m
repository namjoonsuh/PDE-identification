%%%%%%%%%% Generate Parabolic Example Data %%%%%%%%%%%%%%%%
%%%%%%%%%% Roy  Jan.2.2020 %%%%%%%%
%%%% Zero boundary condition
%%%% FTCS scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = ParabolicGen(f,dx,dt,xMax,tMax,nu)
M = floor(xMax/dx);
N = floor(tMax/dt);
U = zeros(N+1,M+1);
U(1,:) = f(0:dx:dx*M);
R = nu*dt/dx^2;
if R > 0.5
    warning(sprintf('Your choice of grid leads to unstable scheme: R=%4f',R))
end
for n = 2:N+1
    U(n,:) = (1-2*R)*U(n-1,:)+...
        R*(U(n-1,[2:end,1])+U(n-1,[end,1:end-1]));
end

end