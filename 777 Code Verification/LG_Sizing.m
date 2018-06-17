function [LG] = LG_Sizing(~,Weight,~,Geom)

%CG WILL BE A FUNCTION, USE VALUES IN INPUT
[~, ~, ~, CG_TO] = Center_of_gravity_new(0, 0, 0, Geom, Weight);
X_cg = CG_TO(1);
Y_cg = CG_TO(2);
Z_cg = CG_TO(3);

aft = 0;    %need to calc, aft shift
fwd = 0;    %need to calc, fwd shift
% B   = 97*12 + 10*12;                        %dist btwn LG
Hcg = (21*12+2) + Y_cg;             %height of CG
nMain = 12; %# main wheels
nNose = 2; %# nose wheels

Z_cg = abs(Z_cg) - (28*12 + 9.5);
Na = abs(Z_cg) + aft;
Nf = abs(Z_cg) - fwd;
% Ma = B - abs(Z_cg) - aft;
% Mf = B - abs(Z_cg) + fwd;
Ma_B = 0.051;
Mf_B = 0.199;

B = (abs(Z_cg) + aft)/(1-Ma_B);

B_f = (abs(Z_cg) - fwd)/(1-Mf_B);
Ma = B - abs(Z_cg) - aft;
Mf = B - abs(Z_cg) + fwd;

MaxSL_main = Weight.MTOW*Na/B;
MaxSL_nose = Weight.MTOW*Mf/B;
MinSL_nose = Weight.MTOW*Ma/B;
dynBrake_nose = 0.31*abs(Hcg)/B*Weight.MTOW;


%% Wheel sizing
Ww_main = MaxSL_main/nMain*1.07;     %wheel load, 1.07 = FAR25 safety factor

% Constants from Raymer, Transport AC
A_Diam = 1.63;
B_Diam = 0.315;
A_Wdth = 0.1043;
B_Wdth = 0.480;

Diam_main = A_Diam*Ww_main^B_Diam;
Wdth_main = A_Wdth*Ww_main^B_Wdth;

LG.Diam_main = Diam_main;
LG.Wdth_main = Wdth_main;
LG.Diam_nose = Diam_main*0.6;   %check this percentage
LG.Wdth_nose = Wdth_main*0.6;   %check this value also gets adjusted


end

