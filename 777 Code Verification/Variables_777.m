function [Req, Area, Main, Geom] = Variables_777(~)
%% Requirements
Req.Pilots        = 2;
Req.Crew          = 12;
Req.Passengers    = 400;
Req.Business      = 49;
Req.Economy       = 351;
Req.Airport_Range = 9395;

%% Area Values
Area.Fuse       = 13125;% Fuselage SWet [ft^2]
Area.Wing       = 4605;
Area.HT         = 1090;% Horizontal Stabilizer Ref. Area [ft^2]
Area.VT         = 573;% Vertical Stabilizer Ref. Area [ft^2]
Area.Wing_CS    = 1.0;% Area of the control surface on wing
Area.Wing_no_CS = Area.Wing*1.0;% Area of no control surface on wing
Area.Nac        = 600;% Nacelle SWet [ft^2]

%% Geometery
% Wing
Geom.Main.sweep    = 33.7;                         %main section sweep [deg]
Geom.Main.dihedral = 6;                            %main section dihedral [deg]
Geom.Main.taper    = 0.114;                         %Taper Ratio
Geom.Main.t_c      = 0.09;                          %Thickness-to-chord ratio
Geom.Main.AR       = 8.68;                           %Aspect Ratio
Geom.Main.Area     = Area.Wing;                    %Area of Wing [ft^2]
Geom.Main.Span     = 212+7/12; %Span [ft]
Geom.Main.x_c      = 0.36;                         %Chordwide location of max thickness

fun                = @(y) (2*Area.Wing/((1+Geom.Main.taper)*Geom.Main.Span)...
                     *(1 - 2*y*(1-Geom.Main.taper)/Geom.Main.Span)).^2;
Geom.Main.MAC      = 2/Area.Wing*integral(@(y)fun(y),0,Geom.Main.Span/2);% Mean Aerodyncamic Chord [ft]

taper                = Geom.Main.taper;
Geom.Main.Chord_Root = Geom.Main.MAC*(3/2)*(1+taper)/(1+taper+taper^2);
Geom.Main.Y_bar      = (Geom.Main.Span/6)*((1+2*taper)/(1+taper));
Geom.Main.Wing_Dist  = 70.7;                   %Wing Tip Location from Nose [ft]
Geom.Tip.sweep       = 40;                        %tip section sweep [deg]
Geom.Tip.dihedral    = 6;                         %tip section dihedral [deg]

% Fuselage
Geom.Fuse.length  = 209+1/12;                     % Fuselage Length [ft] 
Geom.Fuse.Diam    = 20+4/12;                 % Fuselage Diameter [ft]
Geom.Fuse.upsweep = 6*pi/180;               % Fuselage Upsweep [rad]

% Horizontal Stabilizer
Geom.HT.t_c   = 0.09;    %Thickness-to-chord ratio
Geom.HT.x_c   = 0.4;     %Chordwide location of max thickness
Geom.HT.sweep = 35.2;    %Sweep
Geom.HT.AR    = 5.7;     %Aspect Ratio
Geom.HT.Taper = 0.28;    %Taper Ratio
c_Ht          = 1;       %Raymer Constant
L_Ht          = 114.9;
Area.HT       = c_Ht*Geom.Main.MAC*Area.Wing/L_Ht;
Geom.HT.Area  = Area.HT; %Area [ft^2] 
Geom.HT.Span  = sqrt(Geom.HT.AR*Geom.HT.Area); %Wing Area [ft^2]
Geom.HT.Dist  = Geom.Fuse.length - 449/12; %Horizontal Tail (TIP) Location from Nose [ft]

fun = @(y) (2*Geom.HT.Area/((1+Geom.HT.Taper)*Geom.HT.Span)...
           *(1 - 2*y*(1-Geom.HT.Taper)/Geom.HT.Span)).^2;
Geom.HT.MAC = 2/Geom.HT.Area*integral(@(y)fun(y),0,Geom.HT.Span/2);% Mean Aerodyncamic Chord [ft]

taper              = Geom.HT.Taper;
Geom.HT.Chord_Root = Geom.HT.MAC*(3/2)*(1+taper)/(1+taper+taper^2);
Geom.HT.Y_bar      = (Geom.HT.Span/6)*((1+2*taper)/(1+taper));

% Vertical Tail
Geom.VT.t_c   = 0.09;      %Thickness-to-chord ratio
Geom.VT.x_c   = 0.4;       %Chordwide location of max thickness
Geom.VT.sweep = 41.6;      %Sweep
Geom.VT.AR    = 4.0;       %Aspect Ratio
Geom.VT.Taper = 0.29;      %Taper Ratio
c_vt          = 0.09;      %Raymer Constant]
L_vt          = 104.5;
Area.VT       = c_vt*Geom.Main.Span*Area.Wing/L_vt;
Geom.VT.Area  = Area.VT;   %Area [ft^2] 
Geom.VT.Span  = sqrt(Geom.VT.AR*Geom.VT.Area); %Wing Area [ft^2]

fun = @(y) (2*Geom.VT.Area/((1+Geom.VT.Taper)*Geom.VT.Span)...
           *(1 - 2*y*(1-Geom.VT.Taper)/Geom.VT.Span)).^2;
Geom.VT.MAC = 2/Geom.VT.Area*integral(@(y)fun(y),0,Geom.VT.Span/2);% Mean Aerodyncamic Chord [ft]

taper              = Geom.VT.Taper;
Geom.VT.Chord_Root = Geom.VT.MAC*(3/2)*(1+taper)/(1+taper+taper^2);
Geom.VT.Z_bar      = (Geom.VT.Span/6)*((1+2*taper)/(1+taper));

% Nacelle
Geom.Nac.length = 180/12;  % Nacelle "Length" [ft]
Geom.Nac.Diam   = 13;  % Nacelle "Diameter" [ft]

% Landing Gear
Geom.LG.X_pos = Geom.Main.Wing_Dist + Geom.Main.Y_bar*tand(Geom.Main.sweep) + 0.8*Geom.Main.MAC;   % X-Position of Main Landing Gear [ft]
Geom.LG.Z_pos = 230/12;    % Z-Position of Main Landing Gear [ft]
Geom.LG.Z_pos_fuse = 230/12 - 125/12;

%% Main Variables
Main.Range            = 9380;% Maximum Range of Aircraft in [nm]
Main.E                = 0;% Endurance [hr]
Main.V_climb          = 280;% Climb Speed [knots]
Main.Alt              = 35000;% Cruise Altitude [ft]
Main.Climb_rate       = 2500*60;% Climb Rate [fphr]
Main.Thrust           = 115300;
Main.Mach             = 0.84;% Cruise Mach Number
Main.avg_human_weight = 180;% Average Human Weight [lbs]
Main.luggage_weight   = 60;% Average Passenger Luggage Weight [lbs]
Main.S_Land           = 9843;% Maximum landing distance (ft) for (BER) 
Main.rho_rho_SL       = 0.740;% Determined for worst case location pressure ratio (JNB)
Main.Wing_Area        = Area.Wing;% Wing Area [ft^2]
Main.AR               = Geom.Main.AR;
Main.Span             = sqrt(Main.AR*Main.Wing_Area);% Wing Span [ft]
Main.taper            = Geom.Main.taper;% Taper Ratio
Main.sweep            = Geom.Main.sweep;% Sweep
Main.t_c              = Geom.Main.t_c;% Thickness-to-chord ratio
Main.N_z              = 2.5;% Wing Loading
Main.Composite        = 1;% Composite Savings
Main.L_D              = 18;% Lift to Drag Ratio
Main.CL_max           = 2.287;
Main.flight_path      = 5; %Climb flight path angle

[T,P,rho, mu]         = altitude(Main.Alt);
Main.rho_c_rho_SL     = rho/1.225;% Pressure ratio
Main.Mu               = mu; % [lb-s/ft^2]

% Cruise Velocity 
gamma            = 1.4;
R                = 287.1;% [J/Kg-K]
V                = Main.Mach*sqrt(gamma*R*T)*1.9438;%Aircraft Speed [knots]
Main.V_Cruise    = V;  
Main.Thrust_Alt  = Main.Thrust*(Main.rho_c_rho_SL)^0.6;
Main.SFC_Alt     = TSFC_Calc(Main.Alt,Main.Mach);
end