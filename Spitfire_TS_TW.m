function [Constraint, Constraint2] = Spitfire_TS_TW(Main,C_D, C_L)
%% Variables
Constraint = zeros(6,3);
Constraint2 = zeros(1,3);
%% Drag Polar Estimate

% Variables
e_clean = 0.91;% Span Eff. Factor - Clean  
e_TF    = 0.86;% Span Eff. Factor - Takeoff flaps                            
AR      = Main.AR;% Aircraft Aspect Ratio
C_D0_clean    = C_D.clean;% Parasite Drag Coeff. - Clean
C_D0_TF_GUP   = C_D.clean + 0.010;% Parasite Drag Coeff. - Takeoff flaps, gear up
C_D0_TF_GDOWN = C_D.clean + 0.025;% Parasite Drag Coeff. - Takeoff flaps, gear down

%% Landing Field Length
% Variables
S_a = 1000;% Obstacle-clearance distance (ft)
S_Land = Main.S_Land;% Maximum landing distance (ft) for (BER) 
rho_rho_SL = Main.rho_rho_SL;% Determined for worse case location pressure ratio (JNB)
Rho_Cruise = Main.rho_c_rho_SL*1.225*0.00194;%[slug/ft^3]
V_Cruise   = Main.V_Cruise*1.688;%[ft/sec]

q = 1/2*Rho_Cruise*V_Cruise^2;

% Equations
C_L_cruise = C_L;
delta_C_L = 1.56;% Historicl Data
C_L_max = C_L_cruise + delta_C_L;

A = ((rho_rho_SL)/80) * (S_Land - S_a)/1.67 * C_L_max/0.75;
Constraint(1,:) = [1, 0, -A];

%% Takeoff Field Length
% Variables
BFL = S_Land;% Balanced Field Factor

% Equations
TOP = BFL/37.5;% Takeoff Parameter 
B = 1/((rho_rho_SL)*C_L_max*TOP);

Constraint(2,:) = [B, -1, 0];

%% Takeoff Climb
% Variables
G_TF_GDOWN = 0.005;% Flap (DOWN), Gear (Down), Condition (Takeoff)
G_TF_GUP   = 0.024;% Flap (DOWN), Gear (UP), Condition (Takeoff)
G_clean    = 0.012;% Flap (UP), Gear (UP), Condition (Cruise)

% Equations
C = 2*(G_TF_GDOWN + 2*sqrt(C_D0_TF_GDOWN/(pi*AR*e_TF)));%Thrust-to-Weight Ratio - Takeoff flaps down, gear down
Constraint(3,:) = [0, -1, C];

D = 2*(G_TF_GUP + 2*sqrt(C_D0_TF_GUP/(pi*AR*e_TF)));%Thrust-to-Weight Ratio - Takeoff flaps down, gear up
Constraint(4,:) = [0, -1, D];

E = 2*(G_clean + 2*sqrt(C_D0_clean/(pi*AR*e_clean)));%Thrust-to-Weight Ratio - Clean
Constraint(5,:) = [0, -1, E];

%% Ceiling
% Variable

[~,~,rho, ~] = altitude(47000);%Ceiling 47000 [ft]
rho_ceil_rho_SL = rho/1.225;
G_ceiling = 0.001;
c = 0.6;

% Equation
F = 1/(rho_ceil_rho_SL^c)*(G_ceiling + 2*sqrt(C_D0_clean/(pi*AR*e_clean)));
Constraint(6,:) = [0, -1, F];

%% Cruise 
% Variables
rho_c_rho_SL = Main.rho_c_rho_SL;
Constraint2(1,:) = [1/(rho_c_rho_SL^c)*q*C_D0_clean, 1/(rho_c_rho_SL^c)*(1/(q*pi*AR*e_clean)), -1];

end
