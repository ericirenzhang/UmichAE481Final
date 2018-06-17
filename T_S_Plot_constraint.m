function [TW, WS] = T_S_Plot_constraint(Main, Weight,K,W_S, Area, C_L, C_D)
%% Variables
MTOW = Weight.MTOW;% Maximum Takeoff Weight (lbs) 
S    = Main.Wing_Area;% Wing Area [ft^2]

%% Drag Polar Estimate

% Variables
e_clean = 0.91;% Span Eff. Factor - Clean  
e_TF    = 0.86;% Span Eff. Factor - Takeoff flaps                            
e_LF    = 0.81;% Span Eff. Factor - Landing flaps
AR      = Main.AR;% Aircraft Aspect Ratio

C_D0_clean    = C_D.clean;% Parasite Drag Coeff. - Clean
C_D0_TF_GUP   = C_D.clean + 0.020;% Parasite Drag Coeff. - Takeoff flaps, gear up
C_D0_TF_GDOWN = C_D.clean + 0.020 + 0.025;% Parasite Drag Coeff. - Takeoff flaps, gear down

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

if K == 1
    WS = ((rho_rho_SL)/80) * (S_Land - S_a)/1.67 * C_L_max/0.75;
    TW = 0;
end

%% Takeoff Field Length
% Variables
BFL = S_Land;% Balanced Field Factor

% Equations
TOP = BFL/37.5;% Takeoff Parameter 
if K == 2
    TW = W_S./((rho_rho_SL)*C_L_max*TOP);% Thrust-to-Weight Ratio at Takeoff
    WS = 0;
end

%% Takeoff Climb
% Variables
G_TF_GDOWN = 0.005;% Flap (DOWN), Gear (Down), Condition (Takeoff)
G_TF_GUP   = 0.024;% Flap (DOWN), Gear (UP), Condition (Takeoff)
G_clean    = 0.012;% Flap (UP), Gear (UP), Condition (Cruise)

% Equations
if K == 3
    TW = 2*(G_TF_GDOWN + 2*sqrt(C_D0_TF_GDOWN/(pi*AR*e_TF)));%Thrust-to-Weight Ratio - Takeoff flaps down, gear down
    WS = 0;
end

if K == 4
    TW   = 2*(G_TF_GUP + 2*sqrt(C_D0_TF_GUP/(pi*AR*e_TF)));%Thrust-to-Weight Ratio - Takeoff flaps down, gear up
    WS = 0;
end

if K == 5
    TW    = 2*(G_clean + 2*sqrt(C_D0_clean/(pi*AR*e_clean)));%Thrust-to-Weight Ratio - Clean
    WS = 0;
end
          
%% Ceiling
% Variable
rho_c_rho_SL = Main.rho_c_rho_SL;% Pressure ratio at 43,000 ft
G_ceiling = 0.001;

% Equation
if K == 6
    TW= 1/(rho_c_rho_SL^0.6)*(G_ceiling + 2*sqrt(C_D0_clean/(pi*AR*e_clean)));
    WS = 0;
end

%% Cruise 
% Equations
if K == 7
    TW = 1/(rho_c_rho_SL^0.6)*(q*C_D0_clean./(W_S) + (W_S).*(1/(q*pi*AR*e_clean)));
    WS = 0;
end

end
