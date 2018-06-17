function [] = Loads(Main, Area, Weight,figNum)

%%  Initialize Variables
rho_SL = 0.0023769;                 %[slug/ft^3]
sigma = Main.rho_c_rho_SL;
rho_c  = sigma*rho_SL;              %[slug/ft^3]

V_TAS = linspace(0,1000,5000);      %[knots]
V_EAS = linspace(0,500,5000);
CL_cruise_max_pos = 1.5756;    %CHECK THIS!
CL_cruise_max_neg = -0.8301;   %CHECK THIS!!!!!!!!!!!!
CL_cruise = 0.707;

MTOW = Weight.MTOW;                 %[lbs]
S = Area.Wing;                      %[ft^2]

%%   Low Speed
n_LS_pos = linspace(0,2.5,5000);
VEAS_pos = sqrt(n_LS_pos*(2*MTOW/S)/(rho_SL*CL_cruise_max_pos))*0.592483801;     %[knots]

n_LS_neg = linspace(0,-1,5000);
VEAS_neg = sqrt(n_LS_neg*(2*MTOW/S)/(rho_SL*CL_cruise_max_neg))*0.592483801;    %[knots]


%%  Limit Loads
n_Limit_pos = 2.5;
n_Limit_neg = -1;

%%  Stall Line
V_S = sqrt(2*MTOW/(S*rho_SL*CL_cruise_max_pos))*0.592483801;                   %[knots]
V_A = V_S*sqrt(n_Limit_pos);                                                   %[knots]

%%   Gust      CHECK THIS DENSITY and CL_alpha!!!!?!?!?!?!?!?!??!?!?!?!?!??!
b = Main.Span;
Alt = Main.Alt;
AR = Main.AR;
sweep = Main.sweep;
M = Main.Mach;

beta = sqrt(1-M^2);
k = 0.97;
CL_alpha = 2*pi*AR/(2+sqrt(AR^2*(beta/k)^2*(1+tand(sweep)^2/beta^2)+4));
g = 32.2;         %[ft/s^2]

mu = 2*(MTOW/S)/(rho_c*(S/b)*CL_alpha*g);
Kg = 0.88*mu/(5.3+mu);

UB_gust = (66 + (38-66)/30000*(Alt-20000));      %[ft/s]
UC_gust = (50 + (25-50)/30000*(Alt-20000));      %[ft/s]
UD_gust = (25 + (12.5-25)/30000*(Alt-20000));    %[ft/s]

n_GustB_upper = (1 + Kg*CL_alpha*UB_gust*V_EAS*1.688/(498*MTOW/S));
n_GustB_lower = (1 - Kg*CL_alpha*UB_gust*V_EAS*1.688/(498*MTOW/S));
n_GustC_upper = (1 + Kg*CL_alpha*UC_gust*V_EAS*1.688/(498*MTOW/S));
n_GustC_lower = (1 - Kg*CL_alpha*UC_gust*V_EAS*1.688/(498*MTOW/S));
n_GustD_upper = (1 + Kg*CL_alpha*UD_gust*V_EAS*1.688/(498*MTOW/S));
n_GustD_lower = (1 - Kg*CL_alpha*UD_gust*V_EAS*1.688/(498*MTOW/S));

%%  V-n Diagram
figure(figNum);clf;
xlabel('V_{EAS} [knots]','fontsize',16)
ylabel('Load Factor, n','fontsize',16)
hold on

%%   Low Speed
plot(VEAS_pos, n_LS_pos, 'k',VEAS_neg, n_LS_neg, 'k')

%%   Stall
%   Constrain stall n values
[~,i_upper] = min(abs(VEAS_pos-V_S));
n_upper = n_LS_pos(i_upper);
[~,i_lower] = min(abs(VEAS_neg-V_S));
n_lower = n_LS_neg(i_lower);
n_Stall = linspace(n_lower,n_upper,150);
plot(V_S*ones(size(n_Stall)),n_Stall,'k')

%%   Gust Continued

A = rho_SL*CL_cruise_max_pos*S/(2*MTOW)*1.688^2;
B = -Kg*CL_alpha*UB_gust*S/(498*MTOW)*1.688;
C = -1;

VB = (-B + sqrt(B^2 - 4*A*C))/(2*A);

Vel_VB = linspace(0,VB,1000);
n_upper_VB = rho_SL*CL_cruise_max_pos*(Vel_VB*1.68781).^2/(2*MTOW/S);

plot(Vel_VB,n_upper_VB,'k');

[~,i_B_upper] = min(abs(VB - V_EAS));
[~,i_B_lower] = min(abs(VB - V_EAS));
nB_upper = n_GustB_upper(i_B_upper);
nB_lower = n_GustB_lower(i_B_lower);

VC = Main.V_Cruise*sqrt(sigma);         %[knots]
[~,i_C_upper] = min(abs(V_EAS - VC));
[~,i_C_lower] = min(abs(V_EAS - VC));
nC_upper = n_GustC_upper(i_C_upper);
nC_lower = n_GustC_lower(i_C_lower);

VMO = 1.06*VC;                          %[knots]

VD = 1.07*VMO;                         %[knots]
[~,i_D_upper] = min(abs(V_EAS - VD));
[~,i_D_lower] = min(abs(V_EAS - VD));
nD_upper = n_GustD_upper(i_D_upper);
nD_lower = n_GustD_lower(i_D_lower);

plot([VB VC VD VD VC VB],[nB_upper nC_upper nD_upper nD_lower nC_lower nB_lower],'k')


%%   Gust Loads
plot(V_EAS(1:i_B_upper), n_GustB_upper(1:i_B_upper), 'r',...
    V_EAS(1:i_B_lower), n_GustB_lower(1:i_B_lower), 'r',...
    V_EAS(1:i_C_upper), n_GustC_upper(1:i_C_upper), 'c',...
    V_EAS(1:i_C_lower), n_GustC_lower(1:i_C_lower), 'c',...
    V_EAS(1:i_D_upper), n_GustD_upper(1:i_D_upper), 'b',...
    V_EAS(1:i_D_lower), n_GustD_lower(1:i_D_lower), 'b','LineWidth',1)

%%   Limit

V_Limit_pos = linspace(VEAS_pos(end),VD,5000);
V_Limit_neg = linspace(VEAS_neg(end),VC,5000);


plot(V_Limit_pos, n_Limit_pos*ones(size(V_Limit_pos)), 'k',...
    V_Limit_neg, n_Limit_neg*ones(size(V_Limit_neg)), 'k')

plot([VD VD VC],[n_Limit_pos 0 -1],'k')

axis([0 VD+25 -1.5 3.5])



