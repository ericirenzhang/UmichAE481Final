function [Weight, Flight] = Spitfire_Weight(Req, Area, Main, Geom)
%% Variables
% Main Variables
Pilots           = Req.Pilots;
Crew             = Req.Crew;
Passengers       = Req.Passengers;
avg_human_weight = Main.avg_human_weight;% Average Human Weight [lbs]
luggage_weight   = Main.luggage_weight;% Average Passenger Luggage Weight [lbs]
C                = Main.Composite;% Composite Savings
Range            = Main.Range;% Range [nm]
L_D              = Main.L_D;% Lift to Drag Ratio
E                = Main.E;% Endurance [hr]
V_climb          = Main.V_climb;% Climb Speed [knots]
Alt              = Main.Alt;% Cruise Altitude [ft]
Climb_rate       = Main.Climb_rate;% Climb Rate [fphr]
V                = Main.V_Cruise;% Cruise Velocity [knots]
T                = Main.Thrust;% Thrust of MERLIN9X Engine [lbs]

Weight.W_avail_fuel = 0;
Weight.W_trapped    = 0;
Weight.Wf           = 0;
Weight.MTOW         = 0;
Weight.Payload      = 0;
Weight.Crew         = 0;
Weight.We           = 0;
Weight.Wo           = 0;
Weight.Takeoff      = 0;
Weight.Engine       = 0;
Weight.Wing         = 0;
Weight.Fuselage     = 0;
Weight.HT           = 0;
Weight.VT           = 0;
Weight.LD_Gear      = 0;
Weight.Extra        = 0;

% Other Variables
W0 = 700000;%Initial Weight Guess [lbs]
Weight.MTOW = W0; 
tol = 1;
conv = 0;

%% Engine Weight Calculation
W_eng_dry     = 0.521*T^(0.9);% Engine Dry Weight [lbs]
W_eng_oil     = 0.082*T^(0.65);% Engine Oil Weight [lbs]
W_eng_rev     = 0.034*T;% Engine Thrust Reverser Weight [lbs]
W_eng_control = 0.26*T^(0.5);% Engine Control Weight [lbs]
W_eng_start   = 9.33*(W_eng_dry/1000)^(1.078);% Engine Start Weight [lbs]

W_engine = W_eng_dry + W_eng_oil + W_eng_rev + W_eng_control + W_eng_start;% Single Engine Weight [lbs]

%% Fuselage, Wing, Vertical Tail, and Horizontal Tail Weight Calculation
W_Fuse   = 5*Area.Fuse*C;
W_HT     = 5.5*Area.HT*C;
W_VT     = 5.5*Area.VT*C;

%% Weight of Payload and Crew
W_pay  = Passengers*(avg_human_weight + luggage_weight);%Payload Weight [lbs]
W_crew = Pilots*(avg_human_weight + luggage_weight)+...              
         Crew*(avg_human_weight + luggage_weight);%Crew Weight [lbs]
     
%% Fuel Fraction

W1_W0 = 0.970;% Engine start and takeoff
W3_W2 = 0.985;% Climb
W6_W5 = 0.995;% Landing

% Fly to alternate location
Range_alt = 100;%[nm]
SFC_alt   = 0.9; 
V_alt     = 250;%[knots]
L_D_alt   = 12; 
W5_W4     = exp(-Range_alt*SFC_alt/(V_alt*(L_D_alt)));

%% Fuel, Wing, and Empty Weight Calculation
    
i = 1;
while conv == 0    
%% Checks if Weight Fails to converge
    if i > 50
        W0 = 2000000;
        warning('Failed to Converge Weight');
        break
    end

%% Fuel Fraction

% Engine start, warmup and taxi
    Alt_Taxi  = 0; %[ft]
    T_idle    = 8000;%idle thrust [lbf]
    Fuel_idle = TSFC_Calc(Alt_Taxi,0)*T_idle*(15/60);
    
    W1_W0     = (W0-Fuel_idle)/W0; % Engine start, warmup and taxi

% Takeoff
    T_takeoff    = Main.Thrust*2;
    Alt_Takeoff  = 0; %[ft]
    [Temp_Takeoff,~,rho_Takeoff,~] = altitude(Alt_Takeoff);
    e_TF         = 0.86; % Span Eff. Factor - Takeoff flaps                            
    Weight.MTOW  = W1_W0*W0;
    V            = sqrt(Weight.MTOW/(0.5*rho_Takeoff*0.00194*Area.Wing*2.0));
    M_takeoff    = (V/1.9438)/sqrt(1.4*287.1*Temp_Takeoff);
    Fuel_takeoff = TSFC_Calc(Alt_Takeoff,M_takeoff)*T_takeoff*(1/60);

    W2_W1        = (W1_W0*W0-Fuel_takeoff)/(W1_W0*W0);

% Climb and Accel
    Alt_Climb   = 5000; %[ft]
    [Temp_climb,~,~,~] = altitude(Alt_Climb);
    e_clean     = 0.91;% Span Eff. Factor - Clean 
    W           = W2_W1*W1_W0*W0;% [lbs]
    D           = Main.Thrust*2 - W*sind(Main.flight_path);
    Weight.MTOW = W2_W1*W1_W0*W0;
    V_climb     = 540; %[ft/s]
    V_cruise    = Main.V_Cruise*1.9438; %Cruise velocity [ft/s]
    d_he        = (Main.Alt + (V_cruise^2)/2/32.2) - (Alt_Climb + (V_climb^2)/2/32.2);
    M_climb     = (V_climb/1.9438)/sqrt(1.4*287.1*Temp_climb);
    SFC         = TSFC_Calc(Alt_Climb,M_climb);
    
    W3_W2       = exp(-SFC/3600*d_he/(V_climb*(1-D/(Main.Thrust*2))));

    Range_climb  = d_he/sind(Main.flight_path)*0.000164579;%[nm]
    Flight.Range_climb = Range_climb;
    Range_cruise = Range - Range_climb;%[nm]

% Cruise
    Weight.MTOW     = W3_W2*W2_W1*W1_W0*W0;
    [L_D, C_L, C_D] = lift_over_drag(Req, Area, Main, Geom, Weight);
    segsize         = 1000;%Cruise Segment Step size[nm]
    Alt             = Main.Alt;% Altitude in [ft]
    M               = Main.Mach;% Mach at cruise
    W4_W3           = 1;% cruise fuel fraction
    v_cruise        = Main.V_Cruise; %Cruise velocity [knots]
    range_current   = 0;
    
    while range_current < (floor(Range_cruise/(segsize))*segsize)

        Wf_Wi         = exp(-segsize*TSFC_Calc(Alt,M)/(v_cruise*(L_D)));
        W4_W3         = W4_W3*Wf_Wi;
        Wf            = W4_W3*W3_W2*W2_W1*W1_W0*W0;% [lbs]
        P             = 2*(Wf*4.44822162825)/(1.4*(Area.Wing*0.092903)*C_L*M^2);% [N/m^2]
        [Alt]         = Get_Altitude(P);
        [T]           = altitude(Alt);%[K]
        v_cruise      = M*sqrt(1.4*287.1*T)*1.9438;%Aircraft Speed [knots]
        range_current = range_current + segsize;%[nm]

    end

    Range_left      = Range_cruise - floor(Range_cruise/(segsize))*segsize;
    Wf_Wi           = exp(-(Range_left)*TSFC_Calc(Alt,M)/(v_cruise*(L_D)));
    range_current   = range_current + Range_left;%[nm]
    Flight.Increase = Alt - Main.Alt;
    Flight.Ceil     = Alt;
    Flight.Range_cruise = range_current;
    
    W4_W3           = W4_W3*Wf_Wi;

% Final Fuel Fraction
    W6_W0     = W6_W5 * W5_W4 * W4_W3 * W3_W2 * W2_W1 * W1_W0;
    Wf_W0     = (1 - W6_W0)*1.06;% Fuel Fraction

%% Other Weight Calculations
    W_fuel = Wf_W0*W0;%Fuel Weight [lbs]
    W_lg   = 0.043*W0;%Landing Gear Weight [lbs]
    W_xtra = 0.17*W0;%Extra Weight [lbs]

    N_ult  = Main.N_z*1.5; % Ultimate Load Factor
    W_Wing = 0.0051*((W0*N_ult)^0.557)*(Area.Wing^0.649)*...
            (Main.AR^0.5)*(Main.t_c^(-0.4))*((1+Main.taper)^0.1)*...
            secd(Main.sweep)*(Area.Wing_CS^0.1)*C;

%% Weight Derivative
    W_fuel_der = Wf_W0;%Fuel Weight [lbs]
    W_lg_der   = 0.043;%Landing Gear Weight [lbs]
    W_xtra_der = 0.17;%Extra Weight [lbs]
    
    W_Wing_der = 0.0051*(N_ult^(0.557))*(Area.Wing^0.649)*...
            (Main.AR^0.5)*(Main.t_c^(-0.4))*((1+Main.taper)^0.1)*...
            secd(Main.sweep)*(Area.Wing_CS^0.1)*C*(0.557*W0^(0.557-1));
        
    W0_der = W_Wing_der + W_xtra_der + W_lg_der + W_fuel_der;

% Newton Raphson Method
    W0_new = W0 - (2*W_engine + W_Wing + W_HT + W_VT + W_Fuse + W_xtra + W_lg...
        + W_fuel + W_pay + W_crew - W0)/(W0_der - 1);

% Convergence
    if abs(W0_new-W0) <= tol
        conv = 1;
    end
    
    i = i + 1;
    W0 = W0_new;
    Weight.MTOW = W0;
    
end   

%% Results

missionFrac = Req.Airport_Range/Range;

% Real Results
Weight.W_avail_fuel = (1 - W6_W0)*W0*missionFrac;% Available Fuel Weight [lbs]
Weight.W_trapped    = (1 - W6_W0)*0.06*W0;% Weight Trapped Fuel [lbs]
Weight.Wf           = Weight.W_avail_fuel + Weight.W_trapped;% Total Weight of fuel [lbs]
Weight.MTOW         = W0;% Gross Takeoff Weight [lbs]
Weight.Payload      = W_pay;% Weight of Payload (Passenger + luggage) [lbs]
Weight.Crew         = W_crew;% Weight of Crew [lbs]
Weight.We           = W0 - Weight.Wf - Weight.Payload  - Weight.Crew;% Weight Empty [lbs]
Weight.Wo           = Weight.We+Weight.Wf+Weight.Payload+Weight.Crew;% Gross Takeoff Weight [lbs]
Weight.Takeoff      = Weight.We+Weight.Wf+Weight.Payload+Weight.Crew;% Mission Specific Takeoff Weight [lbs]
Weight.Engine       = 2*W_engine;% Weight of 2 Engine [lbs]
Weight.Wing         = W_Wing;% Weight of Wing [lbs]
Weight.Fuselage     = W_Fuse;% Weight of Fuselage [lbs]
Weight.HT           = W_HT;% Weight of Horizontal Tail [lbs]
Weight.VT           = W_VT;% Weight of Vertical Tail [lbs]
Weight.LD_Gear      = W_lg;% Weight of Landing Gear [lbs]
Weight.Extra        = W_xtra;% Weight of All Else Empty [lbs]

end