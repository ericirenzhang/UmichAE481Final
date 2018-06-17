function [T,P,rho, mu] = altitude(Alt)
%% UNITS
% Temperature --> K
% Pressure --> Pa
% density --> kg/m^3
% Dynamic Viscosity --> Pa s
%% Variables
t = 288; %Temperature at sea level (kelvin)
p = 101325; %Pressure at sea level (N/m^2)
R = 287; %Gas constant (J/(kgK))
g = 9.81; %Gravity (m/s^2)
r = 1.225; %(Rho) density at sea level (kg/m^3)
k = 273; %(celcius to kelvin adjustment)
E1 = -6.5/1000; %trophosphere 0-11km (K/km)
E2 = 0; %Tropopause 11-20km (K/km)
E3 = 1/1000; %Stratosphere (1) 20-32km (K/km)

%% Equations

Alt = Alt*0.3048; % Alt in m

if Alt >= 0 && Alt <= 11000
    
    T   = (t + (E1).*(Alt - 0)) - k;%temp for trophosphere
    P   = p.*(((T + k)./t).^(-g/(R*E1)));%pressure for trophosphere
    rho = P/(R.*(T+k));%density for trophosphere
    mu  = (-2.97623003455958e-15*T^4 + -2.24109935584875e-13*T^3 + ...
          7.05814646518648e-12*T^2 + 1.13805307827612e-09*T + ...
          3.56700232977678e-07);
    T   = T + k;
     
elseif Alt > 11000 && Alt <= 20000
    
    T   = (t + E1*11000 + (E2).*(Alt - 11000))- k;%Temp for tropopause
    M   = (p*((t+ E1*11000)/t)^(-g/(R*E1)));%pressure at 11000 m
    P   = M.*(exp(-g.*(Alt - 11000)./(R*(t+E1*11000))));%pressure for tropopause
    rho = P./(R.*(T+k));%density for tropopause
    mu  = (-2.97623003455958e-15*T^4 + -2.24109935584875e-13*T^3 + ...
          7.05814646518648e-12*T^2 + 1.13805307827612e-09*T + ...
          3.56700232977678e-07);
    T   = T + k;
      
elseif Alt > 20000 && Alt <= 30000
    
    T   = (t + E1*11000 + (E3).*(Alt - 20000)) - k;%Temp for stratosphere (1)
    M   = (p*((t+ E1*11000)/t)^(-g/(R*E1)))*(exp(-g*(20000-11000)/(R*(t+E1*11000))));%pressure at 20000 m
    P   = M.*(((T + k)./(t+E1*11000)).^(-g/(R*E3)));%pressure for stratosphere (1)
    rho = P./(R.*(T+k));%density for stratosphere (1)
    mu  = (-2.97623003455958e-15*T^4 + -2.24109935584875e-13*T^3 + ...
          7.05814646518648e-12*T^2 + 1.13805307827612e-09*T + ...
          3.56700232977678e-07);
    T   = T + k; 
    
else

    error('Incorrect Alt: Please enter an altitude between 0 - 98,425 ft')
    
end

end