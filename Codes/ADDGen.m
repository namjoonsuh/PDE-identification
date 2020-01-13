%%%%%%%%%% Generate Advection-Diffusion-D Example Data %%%%%%%%%%%%%%%%
%%%%%%%%%% Roy  Jan.2.2020 %%%%%%%%
%%%% 0 boundary condition
%%%% Lax-Wendroff Diff - ux
%%%% Central Diff - uxx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = ADDGen(f,dx,dt,xMax,tMax,a,nu)
M = floor(xMax/dx);
N = floor(tMax/dt);
U = zeros(N+1,M+1);
U(1,:) = f(0:dx:dx*M);
R1 = a*dt/dx;
R2 = nu*dt/dx^2;
if R1 > 1
    warning(sprintf('Your choice of grid leads to unstable scheme: R1=%4f',R1))
end
for n = 2:N+1
    U(n,2:end-1) = U(n-1,2:end-1)-...
        R1/2*(U(n-1,3:end)-U(n-1,1:end-2))+...
        (R1^2/2+R2)*(U(n-1,3:end)-2*U(n-1,2:end-1)+U(n-1,1:end-2));
        
end
end