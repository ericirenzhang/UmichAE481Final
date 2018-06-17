%% Spitfire Main 
clear all; close all; clc; format longg

%% Optization Initilization
Sguess = 5000;
S_LB   = 4500;
S_UB   = 1e4;

Tguess = 250000;
T_LB   = 1e5;
T_UB   = 4e5;

Mguess = 0.85;
M_LB   = 0.6;
M_UB   = 1;

Altguess = 35000;
Alt_LB   = 20000;
Alt_UB   = 40000;

Sweep = 30;
Sweep_LB = 0;
Sweep_UB = 40;

Taper = 0.2;
Taper_LB = 0.10;
Taper_UB = 0.3;

Wing_Dist = 1000/12;
Wing_LB   = 0;
Wing_UB   = (2039 + 400)/12;

Scale = [5000 200000 0.84 35000 30 0.2 100];

Guess  = [Sguess Tguess Mguess Altguess Sweep Taper Wing_Dist]./Scale;
LB     = [S_LB T_LB M_LB Alt_LB Sweep_LB Taper_LB Wing_LB]./Scale;
UB     = [S_UB T_UB M_UB Alt_UB Sweep_UB Taper_UB Wing_UB]./Scale;

%% Run Optimization FMINCON
options = optimoptions(@fmincon,'PlotFcns',@optimplotfval,'Algorithm','interior-point','Display','iter','DiffMinChange',0);
[X,fval] = fmincon(@ObjectiveFunction,Guess,[],[],[],[],LB,UB,@ConstraintFunction,options);
X = X.*Scale
Obf_fun = fval

%% Get Final Values
% Gets all Variables of Mk35 Aircraft
[Req, Area, Main, Geom] = Variables(X);

% Determines Weight and Flight Parameters
[Weight, Flight] = Spitfire_Weight(Req, Area, Main, Geom, X);

% Determines CG position
[X_CG_TO, X_CG_Empty, x_LG] = Center_of_gravity_new(Geom, Weight);

% Determins Static Margin and Neutral Point
[ SM, xnp, SM_empty ] = Spitfire_Stability(X,X_CG_TO,X_CG_Empty);

% Creates CG Excursion Plot and redifines CG of our aircraft to ensure its
% feasibility
[X_CG_TO, X_CG_Empty, x_LG, CG_TO] = CG_Excursion(Geom, Weight);

% Determines Landing Gear Sizing
[LG] = LG_Sizing(Weight,Geom);

% Determins New Static Margin and Neutral Point
[ SM, xnp, SM_empty ] = Spitfire_Stability(X,X_CG_TO,X_CG_Empty);

%% Plot Final Values
% Creates 3 plots: Trust vs Wing Area, TW-WS, and Drag Polar
Spitfire_Plot( X );

% Creates V-n diagram at cruise Alt
figNum = 777;
Loads(Main, Area, Weight,figNum);
title1 = sprintf('Altitude = %1.0f ft.',X(4));
title(title1,'fontsize',16)

% Creates V-n diagram at 20000 ft (worse case)
figNum = 778;
X(4) = 20000;
[Req, Area, Main, Geom] = Variables(X);
Loads(Main, Area, Weight,figNum);
title1 = sprintf('Altitude = %1.0f ft.',X(4));
title(title1,'fontsize',16)

