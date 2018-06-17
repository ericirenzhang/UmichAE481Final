function [ Cost ] = CostCOC(Weight, Req, Main)
%% Cash Operating Cost (COC)
    % Created Sept 22, 2014
    % Last Revision Oct 12, 2014

%% Direct Operating Cost (Typically: USD/hr or USD/nmi) 

    % ---------------------------------------------------------------------
    % Main Variables
    MTOW = Weight.MTOW;% Maximum Takeoff Weight (lbs)
    
    % ---------------------------------------------------------------------
    % Variables - Pilots
    [ CEF ] = CEF_function( 1995, 2014 );

    % Equations - Pilots
    C_crew = Req.Pilots*(150.0 + 0.0532*(MTOW/1000))*CEF*Main.Airport_t_b;% Crew Cost (USD) - Eqn From Liebeck 1995

    % ---------------------------------------------------------------------
    % Variables - Attendants
    [ CEF ] = CEF_function( 1995, 2014 );

    % Equations - Attendants
    C_attd = 60*Req.Crew*CEF*Main.Airport_t_b;% Cost Flight Attendant (USD)

    % ---------------------------------------------------------------------
    % Variables - Fuel
    P_f = 6;% Fuel Price (USD/gal)
    rho_f = 6.70;% Fuel Density (lbs/gal) Jet A (From Wikipedia)

    % Equations - Fuel
    C_fuel = 1.02*Weight.Wf*(P_f/rho_f);% Fuel Cost (USD)

    % ---------------------------------------------------------------------
    % Variables - Electricity
    W_b = 0;% Weight of Battery (lbs)
    P_elec = 0.47;% Price of Electricity (per/KWh)
    rho_elec = 128;% Energy Density (Wh/kg)

    % Equations - Electricity
    C_elec = 1.05*W_b*P_elec*rho_elec;% Electricity Cost (USD)

    % ---------------------------------------------------------------------
    % Variables - Oil
    P_oil = 62.6;% Oil Price (USD/gal) 
    rho_oil = 8.345;% Oil Density (lbs/gal)

    % Equations - Oil
    W_oil = 0.0125*Weight.Wf*(Main.Airport_t_b/100);% Weight of Oil (lbs)
    C_oil = 1.02*W_oil*(P_oil/rho_oil);% Cost Oil (USD)

    % ---------------------------------------------------------------------
    % Variables - Landing Fees
    [ CEF ] = CEF_function( 1995, 2014 );

    % Equations - Landing Fees
    C_airport = 1.5*((Weight.Takeoff - Weight.Wf)/1000)*CEF;% Landing Fuel Cost (USD)

    % ---------------------------------------------------------------------
    % Variables - Navigation fees
    R = Req.Airport_Range;% Range (nmi)
    [ CEF ] = CEF_function( 1989, 2014 );

    % Equations - Navigation fees                                              
    C_navigation = 0.5*(CEF)*(1.852*R/Main.Airport_t_b)*...
                   sqrt(0.00045359237*MTOW/50);% Navigation Cost (USD) 

    % ---------------------------------------------------------------------
    % Variables - Aircraft maintenance
    R_L = 20;% Maintenance Labor Rate (USD/hr)
    [ CEF ] = CEF_function( 1989, 2014 );

    % Equations - Aircraft maintenance
    C_aircraft = 10.^(3.3191 + 0.8043*log10(MTOW))*CEF;% Total Aircraft Cost (USD)
    C_engines = 10.^(2.3044 + 0.8858.*log10(MTOW))*CEF;% Total Engine Cost (USD)
    C_airframe = C_aircraft - C_engines;% Total Airframe Cost (USD)

    C_ML = 1.03*(3 + 0.067*MTOW/1000)*R_L;% Airframe Labor Cost (USD) 
    C_MM = 1.03*(30*CEF) + 0.79*10^-5 * C_airframe;% Airframe Material Cost (USD)      
    C_airframe_maintenance = (C_ML + C_MM) * Main.Airport_t_b;% Airframe Maintenance Cost (USD)

    % ---------------------------------------------------------------------
    % Variables - Engine Maintenance
    t_bo = 25000;% Time Between Engine Overhauls (hrs)
    T_0 = Main.Thrust;% Engine Maxime Thrust (lbs) 
    n_engines = 2;% Number of Engines

    % Equations - Engine Maintenance
    C_ML = 1.34 * (0.7180 + 0.0317*T_0/1000) * (1100/t_bo + 0.1) * R_L;% Engine Labor Cost (USD)   
    C_MM = 1.34 * (5.43 * 10^-5 * (0.2 * C_engines) - 0.47) ./ ...
           (0.021 * (t_bo/100) + 0.769);% Engine Material Cost (USD)
    C_engine_maintenace = n_engines * (C_ML + C_MM) * Main.Airport_t_b;% Engine Maintenance Cost (USD)

    
    % ---------------------------------------------------------------------
    % Equations - Registration Taxes
    COC = C_crew + C_attd + C_fuel + C_elec + C_oil + C_airport + ...
          C_navigation + C_airframe_maintenance + C_engine_maintenace;% Total Cash Operating Cost
    COC_PassMi = COC/(Req.Passengers*Req.Airport_Range);
    
%% Results

    Cost.COC                  = COC;
    Cost.Crew                 = C_crew;
    Cost.Attendants           = C_attd;
    Cost.Fuel                 = C_fuel;
    Cost.Oil                  = C_oil;
    Cost.Airport_Fees         = C_airport;
    Cost.Navigation_Fees      = C_navigation;
    Cost.Airframe_Maintenance = C_airframe_maintenance;
    Cost.Engine_Maintenance   = C_engine_maintenace;
    Cost.Aircraft             = C_aircraft;
    Cost.Engine               = C_engines;
    Cost.COC_PassMi           = COC_PassMi;

end