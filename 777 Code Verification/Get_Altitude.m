function [Alt] = Get_Altitude(P_guess)

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
E1 = -6.5/1000; %trophosphere 0-11km (K/km)
E2 = 0; %Tropopause 11-20km (K/km)
E3 = 1/1000; %Stratosphere (1) 20-32km (K/km)


%% Equations
if P_guess <= p && P_guess >= 22700
    
    T   = t*(P_guess/p)^(-R*E1/g);
    Alt = (T - t)/E1;%Altitude [m]
    Alt = Alt/0.3048;%Altitude [ft]
     
elseif P_guess < 22700 && P_guess >= 5529
    
    M   = (p*((t+ E1*11000)/t)^(-g/(R*E1)));%pressure at 11000 m
    Alt = -log(P_guess/M)*R*(t+E1*11000)/g + 11000;%Altitude [m]
    Alt = Alt/0.3048;%Altitude [ft]
    
elseif P_guess < 5529 && P_guess >= 1197

    M   = (p*((t+ E1*11000)/t)^(-g/(R*E1)))*(exp(-g*(20000-11000)/(R*(t+E1*11000))));%pressure at 20000 m
    T   = ((P_guess/M)^(-R*E3/g))*(t+E1*11000);
    Alt = (T - t - E1*11000)/E3 + 20000;
    Alt = Alt/0.3048;%Altitude [ft]
    
else
    
    P_guess
    error('Incorrect Pressure: Please enter a pressure between 1197 - 101325 Pa')
        
end

end