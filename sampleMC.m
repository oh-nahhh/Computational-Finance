%Group 1 - Assignment 3
%MC	– Error	as	a	function	of	number	of	paths

S0 = 14; %Initial Stock Price
T = 0.5; %Time to maturity
K = 15; %Strike Price
r = 0.1; %interest rate
sigma = 0.25; %volatility
gamma = 1; %Elasticity
N = [1 10 100 1000 10000 100000]; %nbr of sample paths
M = 10000; %fixed nbr of time steps
nbrSim = length(N); %nbr of simulations
V = zeros(M,1);
exact = bsexact(sigma, r, K, T, S0);
error = zeros(1, nbrSim);
for i = 1:nbrSim
  PartN = N(i);

  deltaT = T/PartN;
  N1 = deltaT:deltaT:T;
  dw = sqrt(deltaT); %Part of the Weiner process

  for n = 1:PartN
    S_old = S0;
    %Eulers method
    for t = N1
      S_new = S_old(:,1) + r*S_old(:,1)*deltaT + sigma*(S_old(:,1)^gamma)*randn*dw;
      S_old = S_new;
    end
    V(n,1) = max(S_new - K, 0);
  end
  Vhat = exp((-r)*T)*sum(V(:,1))/ PartN;
  error(i) = abs(Vhat - exact);
end

h = N ;
Err = error;
h = h(:);
Err = Err(:);
ntest = length(h);
loglog(h,Err,'o-')
title('log-log plot of the sample error for Monte Carlo simulation')

%least square fit
B_q = ones(ntest,2);
B_q(:,2) = log(h);
c_q = log(Err);
L_q = B_q\c_q;
K = L_q(1);
p = L_q(2);

%plot
hold on
error1 = exp(K)*h.^p;
loglog(h,error1,'r')
legend('error points', 'least squares fit','Location','NorthEast')
xlabel('Nbr of Paths');
ylabel('Sample Error');
hold off
