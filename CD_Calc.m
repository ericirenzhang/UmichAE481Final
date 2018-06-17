function [C_D0] = CD_Calc(Area, Main, Geom)

%% Drag Polar
% Variables

e_clean = 0.91;% Span Eff. Factor - Clean 
e_TF    = 0.86;% Span Eff. Factor - Takeoff flaps                            
e_LF    = 0.81;% Span Eff. Factor - Landing flaps

kappa   = 0.95;
t_c     = Main.t_c;
sweep   = Main.sweep;
M       = Main.Mach;

% Cruise Velocity
Rho_Cruise = Main.rho_c_rho_SL*1.225*0.00194;%[slug/ft^3]
S          = Main.Wing_Area;%[ft^2]
V_Cruise   = Main.V_Cruise*1.688;%[ft/sec]

% Parasite Drag
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
    Cf.Wing = 0.1*1.328/sqrt(Re.Wing) + ...
              0.9*0.455/((log10(Re.Wing))^2.58 * (1 + 0.144*M^2)^0.65);
    FF.Wing = (1 + 0.6*Geom.Main.t_c/Geom.Main.x_c + 100*Geom.Main.t_c^4)*...
              (1.32*M^0.18 * cosd(Geom.Main.sweep)^0.28);
    Q.Wing  = 1.1;      

    % Horizontal Tail
    Re.HT   = (Rho_Cruise*V_Cruise*Geom.HT.MAC)/mu;
    Cf.HT   = 0.1*1.328/sqrt(Re.HT) + ...
              0.9*0.455/((log10(Re.HT))^2.58 * (1 + 0.144*M^2)^0.65);
    FF.HT   = (1 + 0.6*Geom.HT.t_c/Geom.HT.x_c + 100*Geom.HT.x_c^4)*...
              (1.32*M^0.18 * cosd(Geom.HT.sweep)^0.28);
    Q.HT    = 1.05;  

    % Vertical Tail
    Re.VT   = (Rho_Cruise*V_Cruise*Geom.VT.MAC)/mu;
    Cf.VT   = 0.1*1.328/sqrt(Re.VT) + ...
              0.9*0.455/((log10(Re.VT))^2.58 * (1 + 0.144*M^2)^0.65);
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
    C_D.mis = (1/Main.Wing_Area)*(D_q.Fuse);

    % Leak and Protuberance Drag
    C_D.LP = 0;

    % Parasitic Drag Calculation
    C_D.clean = (1/Main.Wing_Area)*(Cf.Fuse * FF.Fuse * Q.Fuse * Area.Fuse + ...
                 Cf.Wing * FF.Wing * Q.Wing * (2 * Area.Wing) + ...
                 Cf.HT * FF.HT * Q.HT * (2 * Area.HT) + ...
                 Cf.VT * FF.VT * Q.VT * (2 * Area.VT) + ...
                 2*(Cf.Nac * FF.Nac * Q.Nac * Area.Nac)) + ...
                 C_D.mis + C_D.LP;

% Flap Drag
cf_c = 0.1;
Sf_Sw = 0.9;
df_TO = 10;
df_LDG = 30;

C_D.flap_TO  = 0.9*(cf_c)^1.38*(Sf_Sw)*sind(df_TO).^2;
C_D.flap_LDG = 0.9*(cf_c)^1.38*(Sf_Sw)*sind(df_LDG).^2;

% Landing Gear Drag
D_q_nose = 2*0.25 + 1.0;
D_q_main = 4*(0.25 + 0.15 + 0.15) + 2*(1.4);
C_D.gear = (D_q_nose + D_q_main)/(Main.Wing_Area*0.304);


C_D0.climb    = C_D.clean;
C_D0.TO_Gear  = C_D.clean + C_D.flap_TO + C_D.gear;
C_D0.TO       = C_D.clean + C_D.flap_TO;
C_D0.LDG_Gear = C_D.clean + C_D.flap_LDG + C_D.gear;
C_D0.LDG      = C_D.clean + C_D.flap_LDG;

