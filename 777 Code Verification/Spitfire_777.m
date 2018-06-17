%% Spitfire 777 Verification 
clear all; close all; clc; format longg

% Sets aircraft variables of 777
[Req, Area, Main, Geom] = Variables_777(0);

% Determines Weight and Flight Parameters
[Weight, Flight] = Spitfire_Weight(Req, Area, Main, Geom);

% Determines Landing Gear Sizing
[LG] = LG_Sizing(0,Weight,0,Geom);

% Determines Lift over Drag
[L_D, C_L, C_D] = lift_over_drag(Req, Area, Main, Geom, Weight);

rho_fuel = 6.7; %[lbs/gal]
Fuel_Capacity = Weight.W_avail_fuel/rho_fuel;%[gal]