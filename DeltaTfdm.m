%FDM	- Error	as	function	of	time	step
%Group 1 - Assignment 3

clear all; close all;

%Parameters
K = 15; % Strike price
r = 0.1; %Interest rate
sigma = 0.25; % volatility
T = 0.5; % Maturity time
s_max = 4*K; % Max stock price
M = 1000; % nbr of price steps
%N = 10000; % nbr of time steps
nbrN = 20:10:250; %convergence points
sim = length(nbrN);
for i = 1:sim
  [v,S] = implicitFDM(K,r,sigma,T,M,nbrN(i));
  for j = 1:M-1
      v_exact = bsexact(sigma, r, K, T, S(j));
      Error(j) = abs(v(j,1)-v_exact); % error
  end
  MeanError(i) = mean(Error);
end

h = nbrN;
err = MeanError;
h = h(:);
err = err(:);
ntest = length(h);
loglog(h,err,'o-')
title('log-log plot of varying the time step for FDM')

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
Error1= exp(K)*h.^p;
loglog(h,Error1,'r')
legend('errors', 'least squares fit','Location','NorthEast')
xlabel('Nbr of time points');
ylabel('Error');
hold off
