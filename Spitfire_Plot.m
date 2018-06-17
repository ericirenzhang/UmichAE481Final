function [Req, Area, Main ] = Spitfire_Plot( X )
%% T-S Plot with COC Contour Overlay
[Req, Area, Main, Geom] = Variables(X);
[ObjFun] = T_S_Plot(Req, Area, Main, Geom, X);

figure(2)
hold on
plot(X(1), X(2),'ko','linewidth',5)

%% T-W vs T-S Plot
[Weight, Flight] = Spitfire_Weight(Req, Area, Main, Geom, X);
[L_D, C_L, C_D] = lift_over_drag(Area, Main, Geom, Weight);
Main.L_D = L_D;
[C, C2] = Spitfire_TS_TW(Main,C_D, C_L);

W_S = linspace(0,250,1000);
T_W = linspace(0,0.5,1000);

e_clean = 0.91;% Span Eff. Factor - Clean  
e_TF    = 0.86;% Span Eff. Factor - Takeoff flaps 
AR      = Main.AR;% Aircraft Aspect Ratio
C_D0_clean    = C_D.clean;% Parasite Drag Coeff. - Clean
C_D0_TF_GUP   = C_D.clean + 0.010;% Parasite Drag Coeff. - Takeoff flaps, gear up
C_D0_TF_GDOWN = C_D.clean + 0.025;% Parasite Drag Coeff. - Takeoff flaps, gear down
Rho_Cruise = Main.rho_c_rho_SL*1.225*0.00194;%[slug/ft^3]
V_Cruise   = Main.V_Cruise*1.688;%[ft/sec]
q = 1/2*Rho_Cruise*V_Cruise^2;
G_TF_GDOWN = 0.005;% Flap (DOWN), Gear (Down), Condition (Takeoff)
G_TF_GUP   = 0.024;% Flap (DOWN), Gear (UP), Condition (Takeoff)
G_clean    = 0.012;% Flap (UP), Gear (UP), Condition (Cruise)

T_W1 = linspace(2*(sqrt(4*C_D0_TF_GDOWN/(pi*AR*e_TF)) + G_TF_GDOWN),0.5,1000);
T_W2 = linspace(2*(sqrt(4*C_D0_TF_GUP/(pi*AR*e_TF)) + G_TF_GUP),0.5,1000);
T_W3 = linspace(2*(sqrt(4*C_D0_clean/(pi*AR*e_clean)) + G_clean),0.5,1000);
W_S1 = ((T_W/2 - G_TF_GDOWN) + sqrt((T_W/2 - G_TF_GDOWN).^2 + ...
       (4*C_D0_TF_GDOWN/(pi*AR*e_TF))))./(2/(q*pi*AR*e_TF));
W_S2 = ((T_W/2 - G_TF_GUP) + sqrt((T_W/2 - G_TF_GUP).^2 + ...
       (4*C_D0_TF_GUP/(pi*AR*e_TF))))./(2/(q*pi*AR*e_TF));
W_S3 = ((T_W/2 - G_clean) + sqrt((T_W/2 - G_clean).^2 + ...
       (4*C_D0_clean/(pi*AR*e_clean))))./(2/(q*pi*AR*e_clean));
  

figure()
set(0,'DefaultAxesColorOrder',[0 0 1; 0 0.5 0; 0.75 0.75 0; 1 0 0;...
    0 0 0; 1 0 1; 0 1 1; 0.75 0.75 0; 1 0 0; 0 0 0]);
set(0,'DefaultLineLinewidth',2);

plot(-C(1,3)*ones(1,size(W_S,2)),T_W, W_S,W_S.*C(2,1), W_S,C(3,3)*ones(1,size(T_W,2)),...
     W_S,C(4,3)*ones(1,size(T_W,2)), W_S,C(5,3)*ones(1,size(T_W,2)),...
     W_S,(C(6,3)*ones(1,size(T_W,2))), W_S,(C2(1,1).*(1./W_S) + C2(1,2).*W_S),...
     W_S1,T_W1, W_S2,T_W2, W_S3,T_W3)
hold on
plot(Weight.MTOW/X(1), X(2)/Weight.MTOW,'ko','linewidth',5)

set(gca,'FontSize',16)
legend('Landing','Takeoff','1st Climb','2nd Climb','3rd Climb',...
       'Ceiling','Cruise','Location','EastOutside')
axis([0,250,0,0.5])
xlabel('W/S')
ylabel('T/W')

%% Drag Polar
% Variables
C_L = linspace(0,2,100);
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

% AVL && Induced Drag
% [AVL] = AVLrun(X,W);
% C_L = AVL(2);
% e_clean = AVL(4);
% e_TF    = AVL(4) - 0.05;
% e_LF    = AVL(4) - 0.10

% Wave Drag
M_DD = kappa/cosd(sweep) - t_c/(cosd(sweep))^2 - C_L./(10*(cosd(sweep))^3);
M_crit = M_DD - (0.1/80)^(1/3);

if M > M_crit
    C_D_wave = 20.*(M-M_crit).^4;
else
    C_D_wave = 0;
end

% Total Drag and Drag Polar
% C_D = C_D.clean + C_D_induced + C_D_wave;
C_L0 = 0.0834;
C_D.cruise   = C_D.clean + ((C_L-C_L0).^2)./(pi*Main.AR*e_clean) + C_D_wave;
C_D.TO_Gear  = C_D.clean + ((C_L-C_L0).^2)./(pi*Main.AR*e_TF) + C_D.flap_TO + C_D.gear;
C_D.TO       = C_D.clean + ((C_L-C_L0).^2)./(pi*Main.AR*e_TF) + C_D.flap_TO;
C_D.LDG_Gear = C_D.clean + ((C_L-C_L0).^2)./(pi*Main.AR*e_LF) + C_D.flap_LDG + C_D.gear;
C_D.LDG      = C_D.clean + ((C_L-C_L0).^2)./(pi*Main.AR*e_LF) + C_D.flap_LDG;

L_D = C_L./C_D.cruise;

figure()
subplot(2,1,1)
plot(C_D.cruise,C_L,'k-','linewidth',1)
hold on
plot(C_D.TO_Gear,C_L,'b-','linewidth',1)
hold on
plot(C_D.TO,C_L,'r-','linewidth',1)
hold on
plot(C_D.LDG_Gear,C_L,'g-','linewidth',1)
hold on
plot(C_D.LDG ,C_L,'m-','linewidth',1)
xlabel('C_D')
ylabel('C_L')
title('Drag Polar')
legend('Clean','Takeoff flaps + gear','Takeoff flaps',...
       'Landing flaps + gear','Landing flaps','Location','SouthEast');

subplot(2,1,2)
plot(C_L,L_D,'linewidth',2)
xlabel('C_L')
ylabel('L/D')
end

