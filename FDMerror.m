%Group 1 - Assignment 3

clear all; close all;
%Parameters
K = 15; % Strike price
r = 0.1; %Interest rate
sigma = 0.25; % volatility
T = 0.5; % Maturity time
s_max = 4*K; % Max stock price
%M = 1000; % nbr of price steps
N = 10000; % nbr of time steps

nbrM = 10:10:150; %convergence points
sim = length(nbrM);
for i = 1:sim
  [v,S] = implicitFDM(K,r,sigma,T,nbrM(i),N);
  for j = 1:nbrM(i)-1
    VE = bsexact(sigma, r, K, T, S(j));
    Error(j) = abs(v(j,1)-VE); % error
  end
  Mean_err(i) = mean(Error);
end

h = nbrM;
err = Mean_err;

h = h(:);
err= err(:);
ntest = length(h);
loglog(h,err,'o-')
axis([9 200 0 0.1])
title('log-log plot of varying the price steps for FDM')

%least square fit
B_q = ones(ntest,2);
B_q(:,2) = log(h);
c_q = log(err);
L_q = B_q\c_q;
K = L_q(1);
p = L_q(2);
disp(' ')
disp(fprintf('Least squares fit gives E(h) = %g * h^%g',exp(K),p))
disp(' ')

%plots
hold on
Error1 = exp(K)*h.^p;
loglog(h,Error1,'r')
legend('errors', 'least squares fit','Location','NorthEast')
xlabel('Nbr of price points');
ylabel('Error');

hold off
