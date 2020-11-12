function [v,sv] = implicitFDM(K,r,sigma,T,M,N)
%Allocate
B = zeros(M-1,1);
A = zeros(M-1,M-1);
v = zeros(M-1,N);
gamma = 1.0;
% Euler backward method to price European call
Smax = 4*K; % max share price
ds = Smax/M; % Price step
dt = T/N; % Time step
s = ds:ds:Smax;
t = dt:dt:T;
ss = sigma*sigma; %constant
sv = s;
for i =1:M-1
 v(i, N) = max(s(i)-K,0.0); % Final condition V(S,T)=max(S-K,0)
end
%Perform Euler backwards
n = N;
while n >1
 %Euler backwards in time
 for j = 2:(M-1)
 alpha(j) = - 0.5*ss*s(j)^(2*gamma)*dt/(ds*ds) + 0.5*r*s(j)*dt/ds;
 beta(j) = 1 + ss*s(j)^(2*gamma)*dt/(ds*ds) + dt*r;
 zeta(j) = - 0.5*r*s(j)*dt/ds - 0.5*ss*s(j)^(2*gamma)*dt/(ds*ds);
 end
 %tridiagonal matrix, A
 A(1,1) = 1; %first element of A
 for j = 2:M-2
 A(j,j-1:j+1) = [alpha(j) beta(j) zeta(j)]; %fill the diagonals of
A
 end

A(M-1,M-2:M-1) = [alpha(M-1) beta(M-1)]; %last elements of A
% Final Boundary condition V(S,t)=Smax-Ke^(-r(T-t))
B(M-1) = zeta(M-1) * (Smax - K*exp(-r*(T-t(n-1))));
%Implicit method: v^(n-1) = A*v^(n) + b.
v(:,n-1) = A\(v(:,n) - B);
n = n-1; %count backwards
end
end
