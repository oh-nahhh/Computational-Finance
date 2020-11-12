%Group 1 - Assignment 3
%Monte	Carlo	Method	– used	for	calculating	MC Delta

function [v] = mc(T,K,r,sigma,x,y, z)
%x - nbr of simulation steps
%y = time steps
%z - nbr of sample paths

Smax = 4*K;
ds = Smax/x;
deltaT = T/x;

%initial value of S_y for each simulation step
for j = 1:x
  S_y(:,j) = ds*j;
end
%run for x simulation steps
for j = 1:x
  %S_plus
  for k = 1:z
    s_old(1,1)=S_y(j);
    for i = 1:y
        Z = randn;
        s_new = s_old(:,1) + r.*s_old(:,1).*deltaT + sigma.*s_old(:,1).*Z*sqrt(deltaT);
        s_old= s_new;
    end
    Vp(k,1) = max(s_new(:,1) - K, 0);
  end
  Vp_hat = (exp(-r.*T)).*sum(Vp(:,1))./z;

  %S_xinus
  for k = 1:z
      s_old(1,1)=S_y(j);
      for i = 1:y
        Z = randn;
        s_new = s_old(:,1) + r.*s_old(:,1).*deltaT + sigma.*s_old(:,1).*((-1)*Z).*sqrt(deltaT);
        s_old= s_new;
      end
      Vm(k,1) = max(s_new(:,1) - K, 0);
  end

Vm_hat = (exp(-r.*T)).*sum(Vm(:,1))./z;

%Antithetic xC (V_plus + V_minus )/2
Vhat(j) = mean(Vp_hat + Vm_hat )/2; %Save V_hat for each simulation step
end
v = Vhat;
end
