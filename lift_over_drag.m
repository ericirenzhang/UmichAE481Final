function [L_D, C_L, C_D] = lift_over_drag(Area, Main, Geom, Weight)
%% Variables
W = Weight.MTOW;
L = W;
e_clean = 0.91;% Span Eff. Factor - Clean 
kappa   = 0.95;
t_c     = Main.t_c;
sweep   = Main.sweep;
M       = Main.Mach;

%% Cruise Velocity
Rho_Cruise = Main.rho_c_rho_SL*1.225*0.00194;%[slug/ft^3]
S          = Main.Wing_Area;%[ft^2]
V_Cruise   = Main.V_Cruise*1.688;%[ft/sec]

%% Parasite Drag
mu = Main.Mu; % [lb/(ft-s)]

% Fuselage
Re.Fuse = (Rho_Cruise*V_Cruise*Geom.Fuse.length)/mu;
Cf.Fuse = 0.05*1.328/sqrt(Re.Fuse) + ...
          0.95*0.455/((log10(Re.Fuse))^2.58 * (1 + 0.144*M^2)^0.65);
f       = Geom.Fuse.length/(Geom.Fuse.Diam);
FF.Fuse = (1+ (60/f^3) + f/400); % Form Factor
Q.Fuse  = 1.0; % Interference Factor

% Wing
Re.Wing = (Rho_Cruise*V_Cruise*Geom.Main.MAC)/mu;
Cf.Wing = 0.0*1.328/sqrt(Re.Wing) + ...
          1.0*0.455/((log10(Re.Wing))^2.58 * (1 + 0.144*M^2)^0.65);
FF.Wing = (1 + 0.6*Geom.Main.t_c/Geom.Main.x_c + 100*Geom.Main.t_c^4)*...
          (1.32*M^0.18 * cosd(Geom.Main.sweep)^0.28);
Q.Wing  = 1.1;      

% Horizontal Tail
Re.HT   = (Rho_Cruise*V_Cruise*Geom.HT.MAC)/mu;
Cf.HT   = 0.0*1.328/sqrt(Re.HT) + ...
          1.0*0.455/((log10(Re.HT))^2.58 * (1 + 0.144*M^2)^0.65);
FF.HT   = (1 + 0.6*Geom.HT.t_c/Geom.HT.x_c + 100*Geom.HT.x_c^4)*...
          (1.32*M^0.18 * cosd(Geom.HT.sweep)^0.28);
Q.HT    = 1.05;  

% Vertical Tail
Re.VT   = (Rho_Cruise*V_Cruise*Geom.VT.MAC)/mu;
Cf.VT   = 0.0*1.328/sqrt(Re.VT) + ...
          1.0*0.455/((log10(Re.VT))^2.58 * (1 + 0.144*M^2)^0.65);
FF.VT   = (1 + 0.6*Geom.VT.t_c/Geom.VT.x_c + 100*Geom.VT.x_c^4)*...
          (1.32*M^0.18 * cosd(Geom.VT.sweep)^0.28);
Q.VT    = 1.05;  

% Nacelles
Re.Nac = (Rho_Cruise*V_Cruise*Geom.Nac.length)/mu;
Cf.Nac = 0.455/((log10(Re.Nac))^2.58 * (1 + 0.144*M^2)^0.65);
FF.Nac = 1 + 0.35/(Geom.Nac.length/Geom.Nac.Diam);
Q.Nac  = 1.3;

% Missing Form Drag
D_q.Fuse = 3.83*(Geom.Fuse.upsweep^2.5)*(pi*(Geom.Fuse.Diam/2)^2);
C_D_mis = (1/Main.Wing_Area)*(D_q.Fuse);

% Leak and Protuberance Drag
C_D_LP = 0;

% Parasitic Drag
C_D.clean = (1/Main.Wing_Area)*(Cf.Fuse * FF.Fuse * Q.Fuse * Area.Fuse + ...
             Cf.Wing * FF.Wing * Q.Wing * (2 * Area.Wing) + ...
             Cf.HT * FF.HT * Q.HT * (2 * Area.HT) + ...
             Cf.VT * FF.VT * Q.VT * (2 * Area.VT) + ...
             2*(Cf.Nac * FF.Nac * Q.Nac * Area.Nac)) + ...
             C_D_mis + C_D_LP;

% S_wet = Area.Fuse + 2*(Area.Wing + Area.HT + Area.VT);
% f = C_f*S_wet;% Equivalent Parasite Area (ft^2)
% 
% C_D0_clean    = f/S;% Parasite Drag Coeff. - Clean

%% AVL && Induced Drag
C_L = L*2/(Rho_Cruise*V_Cruise^2*S);

% [AVL] = AVLrun(X,W,C_D.clean,C_L);
% C_L = AVL(2);
% C_D.AVL = AVL(3);
% e_clean = AVL(4);

%% Wave Drag
M_DD = kappa/cosd(sweep) - t_c/(cosd(sweep))^2 - C_L/(10*(cosd(sweep))^3);
M_crit = M_DD - (0.1/80)^(1/3);

if M > M_crit
    C_D_wave = 20*(M-M_crit)^4;
else
    C_D_wave = 0;
end

%% Total Drag and Drag Polar
% C_D = C_D.clean + C_D_induced + C_D_wave;
C_L0 = 0.0834;
C_D.total = C_D.clean + ((C_L-C_L0)^2)/(pi*Main.AR*e_clean) + C_D_wave;
% L_D = C_L/C_D.AVL;
 L_D = C_L/C_D.total;
end