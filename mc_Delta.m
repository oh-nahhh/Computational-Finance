%Group 1 - Assignment 3
%MC	vs.	Exact	Delta

clear all; close all; clc;

M = 100; %nbr of time steps
N =100000;%nbr of sample paths
Steps= 100; %nbr of simulation steps
S0 = 14; %Initial Stock Price
K = 15; %Strike Price
sigma = 0.25; %volatility
r = 0.1; %interest rate
T = 0.5; %Time to maturity
Smax = 4*K;
ds = Smax/M;
S = ds:ds:Smax;
deltaT = T/M; %time step
Vhat = mc(T,K,r,sigma,M,Steps, N);

%calculate Delta by taking the derivative of V_hat at each point
for k = 2:M-1
  Delta(k) = (Vhat(k+1)-Vhat(k-1))/(2*ds); %centered differences
end

%exact solution
d1 = ((log(S./K) + (r+0.5.*sigma.^2).*T)./(sigma.*sqrt(T)))./sqrt(2);
exact = 0.5.*(1+erf(d1));

%Plot against the Stock
plot(S(1:end-1),Delta,'r')
title('MC vs Exact Delta');
xlabel('stock')
hold on
plot(S,exact);
