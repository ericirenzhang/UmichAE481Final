%% Calculates the Dynamic V
T = [-40; -20; 0; 10; 20; 30; 40; 50; 60; 70; 80; 90];
T = (T - 32)./1.80;
mu = [3.29; 3.34; 3.38; 3.44; 3.50; 3.58; 3.60; 3.68; 3.75; 3.82; 3.86; 3.90]*10^-7;
  
P = polyfit(T,mu,4)
mu_new = polyval(P,T);
plot(T,mu,'g--',T,mu_new,'r--')
