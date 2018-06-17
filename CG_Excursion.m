function [X_CG_TO, X_CG_Empty, x_LG, CG_TO] = CG_Excursion(Geom, Weight)
%% CENTER OF GRAVITY
Weightplot = Weight;
%% Weight & Moment of Wing
W_wing = Weight.Wing;% Wing Weight [lbs]

% Wing_Dist = 0:1:Geom.Fuse.length;% Vector of Wing Distance [inch]
Wing_Dist = Geom.Main.Wing_Dist;
x_mac  = Wing_Dist + Geom.Main.Y_bar*tand(Geom.Main.sweep);% Mean Aerodynamic Chord from Nose [ft]
x_Wing = x_mac + 0.25*Geom.Main.MAC;% Wing X-CG Position [ft]
y_Wing = 0;% Wing Y-CG Position [ft]
z_Wing = (43.57 - 70)/12;% Wing Z-CG Position [ft] based on CAD    

%% Weight & Moment of Fuel
W_fuel = Weight.Wf;% Weight of Fuel [lbs]

x_fuel = (x_mac + 0.25*Geom.Main.MAC);% Fuel X-CG Position [ft]
y_fuel = 0;% Fuel Y-CG Position [ft] 
z_fuel = (43.57 - 70)/12;% Fuel Z-CG Position [ft] based on CAD

%% Weight % Moment of Payload
Weight_payload = Weight.Payload;% Weight of Payload [lbs]

x_payload = Geom.Fuse.length*0.5;% Payload CG X-CG Position [ft]
y_payload = 0;% Payload CG Y-CG Position [ft]
z_payload = 0;% Payload CG Z-CG Position [ft]

%% Weight & Moment of Crew
Weight_crew = Weight.Crew;% Weight of Crew [lbs]

x_crew = Geom.Fuse.length*0.5;% Crew CG X-CG Position [ft]
y_crew = 0;% Crew CG Y-CG Position [ft]
z_crew = 0;% Crew CG Z-CG Position [ft]

%% Weight & Moment of Fuselage
W_Fuselage = Weight.Fuselage;% Weight of Fuselage [lbs]

x_Fuselage = Geom.Fuse.length*0.47;% Fuselage CG X-CG Position [ft]
y_Fuselage = 0;% Fuselage CG Y-CG Position [ft]
z_Fuselage = 0;% Fuselage CG Z-CG Position [ft]

%% Weight & Moment of Horizontal Tail
W_HT = Weight.HT;% Weight of Horizontal Tail [lbs]

x_Horizontal_Tail = (2500 + 400)/12;% Horizontal Tail CG X-CG Position [ft] based on CAD
y_Horizontal_Tail = 0;% Horizontal Tail CG Y-CG Position [ft]
z_Horizontal_Tail = (17.72 + 50)/12;% Horizontal Tail CG Z-CG Position [ft] based on CAD

%% Weight & Moment of Vertical Tail
W_VT    = Weight.VT;% Weight of Vertical Tail [lbs]

x_Vertical_Tail = (2400 + 400)/12;% Vertical Tail CG X-CG Position [ft] based on CAD
y_Vertical_Tail = 0;% Vertical Tail CG Y-CG Position [ft] based on CAD
z_Vertical_Tail = (179.29 + 80)/12;% Vertical Tail CG Z-CG Position [ft] based on CAD

%% Weight & Moment of All else Empty
Weight_Extra = Weight.Extra;% Weight of All Else Empty [lbs]

x_All_else_empty = Geom.Fuse.length*0.5;% All Else Empty CG X-CG Position [ft]
y_All_else_empty = 0;% All Else Empty CG Y-CG Position [ft]
z_All_else_empty = 0;% All Else Empty CG Z-CG Position [ft]

%% Weight & Moment of Engine
Weight_Engine = Weight.Engine;% Weight of 2 Engine [lbs]

x_engine = x_mac - (275/2/12);% Engine Empty CG X-CG Position [ft]
y_engine = 0;% Engine Empty CG Y-CG Position [ft]
z_engine = -110/12;% Engine Empty CG Z-CG Position [ft] based on CAD

%% Weight & Moment of Landing Gear
Weight_LG = Weight.LD_Gear;

x_LG = ((Geom.LG.X_pos)*12 + (225/12)*2)/14;% Landing Gear CG X-CG Position [ft]
y_LG = 0;% Landing Gear CG Y-CG Position [ft]
z_LG = -Geom.LG.Z_pos;% Landing Gear CG Z-CG Position [ft]

%% Center of Gravity Position
W_empty = Weight.We;

%% X
X_moment_empty = x_Wing*(W_wing) +...
                 x_Horizontal_Tail*(W_HT) +...
                 x_Vertical_Tail*(W_VT) + ...
                 x_Fuselage*(W_Fuselage) + ...
                 x_engine*(Weight_Engine) + ...
                 x_All_else_empty*(Weight_Extra) + ...
                 x_LG*(Weight_LG);
X_CG_Empty     = X_moment_empty./W_empty;
   
%% Y
Y_moment_empty = y_Wing*(W_wing) +...
                 y_Horizontal_Tail*(W_HT) +...
                 y_Vertical_Tail*(W_VT) + ...
                 y_Fuselage*(W_Fuselage) + ...
                 y_engine*(Weight_Engine) + ...
                 y_All_else_empty*(Weight_Extra) + ...
                 y_LG*(Weight_LG);
Y_CG_Empty     = Y_moment_empty./W_empty;

%% Z
Z_moment_empty = z_Wing*(W_wing) +...
                 z_Horizontal_Tail*(W_HT) +...
                 z_Vertical_Tail*(W_VT) + ...
                 z_Fuselage*(W_Fuselage) + ...
                 z_engine*(Weight_Engine) + ...
                 z_All_else_empty*(Weight_Extra) + ...
                 z_LG*(Weight_LG);
Z_CG_Empty     = Z_moment_empty./W_empty;             
   
%% CG

X_moment_full = X_moment_empty + ...
                x_fuel*W_fuel + ...
                x_payload*Weight_payload + ...
                x_crew*Weight_crew;
X_CG_TO       = X_moment_full./Weight.MTOW;

Y_moment_full = Y_moment_empty + ...
                y_fuel*W_fuel + ...
                y_payload*Weight_payload + ...
                y_crew*Weight_crew;
Y_CG_TO       = Y_moment_full./Weight.MTOW;

Z_moment_full = Z_moment_empty + ...
                z_fuel*W_fuel + ...
                z_payload*Weight_payload + ...
                z_crew*Weight_crew;
Z_CG_TO       = Z_moment_full./Weight.MTOW;

CG_TO = [X_CG_TO, Y_CG_TO, Z_CG_TO];

% Variables
CG = zeros(6,4,3);
Weight = zeros(6,4);
Moment = zeros(6,4,3);

Weight(:,1) = W_empty; % Empty Weight CG
Moment(:,1,1) = X_moment_empty;
Moment(:,1,2) = Y_moment_empty;
Moment(:,1,3) = Z_moment_empty;

moment_comp = [x_fuel*W_fuel, y_fuel*W_fuel, z_fuel*W_fuel;
               x_payload*Weight_payload, y_payload*Weight_payload, z_payload*Weight_payload;
               x_crew*Weight_crew, y_crew*Weight_crew, z_crew*Weight_crew];

Weight_comp = [W_fuel; Weight_payload; Weight_crew];

for i = 1:3
    Moment(i,2,1) = Moment(i,1,1) + moment_comp(i,1);
    Moment(i,2,2) = Moment(i,1,2) + moment_comp(i,2);
    Moment(i,2,3) = Moment(i,1,3) + moment_comp(i,3);
    
    Moment(i+3,2,1) = Moment(i+3,1,1) + moment_comp(i,1);
    Moment(i+3,2,2) = Moment(i+3,1,2) + moment_comp(i,2);
    Moment(i+3,2,3) = Moment(i+3,1,3) + moment_comp(i,3);
    
    Weight(i,2) = Weight(i,1) + Weight_comp(i,1);
    Weight(i+3,2) = Weight(i+3,1) + Weight_comp(i,1);
    
    True = 0;
    
    for j = 1:3
        
        if i ~= j && True == 0;
            Moment(i,3,1) = Moment(i,2,1) + moment_comp(j,1);
            Moment(i,3,2) = Moment(i,2,2) + moment_comp(j,2);
            Moment(i,3,3) = Moment(i,2,3) + moment_comp(j,3);

            Weight(i,3) = Weight(i,2) + Weight_comp(j,1);
            
            for k = 1:3
                if i ~= j && j ~= k && i ~= k 
                    Moment(i,4,1) = Moment(i,3,1) + moment_comp(k,1);
                    Moment(i,4,2) = Moment(i,3,2) + moment_comp(k,2);
                    Moment(i,4,3) = Moment(i,3,3) + moment_comp(k,3);

                    Weight(i,4) = Weight(i,3) + Weight_comp(k,1);
                end
            end
            
            True = 1;
            
        elseif i ~= j && True == 1;
            Moment(i+3,3,1) = Moment(i+3,2,1) + moment_comp(j,1);
            Moment(i+3,3,2) = Moment(i+3,2,2) + moment_comp(j,2);
            Moment(i+3,3,3) = Moment(i+3,2,3) + moment_comp(j,3);

            Weight(i+3,3) = Weight(i+3,2) + Weight_comp(j,1);
            
            for k = 1:3
                if i ~= j && j ~= k && i ~= k 
                    Moment(i+3,4,1) = Moment(i+3,3,1) + moment_comp(k,1);
                    Moment(i+3,4,2) = Moment(i+3,3,2) + moment_comp(k,2);
                    Moment(i+3,4,3) = Moment(i+3,3,3) + moment_comp(k,3);

                    Weight(i+3,4) = Weight(i+3,3) + Weight_comp(k,1);
                end
            end
            
        end
        
    end
    
end

CG(:,:,1) = Moment(:,:,1)./Weight;
CG(:,:,2) = Moment(:,:,2)./Weight;
CG(:,:,3) = Moment(:,:,3)./Weight;
Weight_CG = Weight;

%% Plot CG

Z_cg = CG(end,end,1);
aft = 130.4418;  
fwd = 130.1807; 
Weight_min = min(min(abs(Weight_CG(:,:))));
Weight_max = max(max(abs(Weight_CG(:,:))));
Weight_plot = linspace(Weight_min-1e5,Weight_max+1e5,100);

figure(22)
plot(130.3132,Weightplot.We,'ro','linewidth',4)
hold on
plot(130.2768,Weightplot.Wo,'go','linewidth',4)
hold on
plot(CG(1,:,1),Weight_CG(1,:),'k')
hold on
plot(CG(2,:,1),Weight_CG(2,:),'k')
hold on
plot(CG(3,:,1),Weight_CG(3,:),'k')
hold on
plot(CG(4,:,1),Weight_CG(4,:),'k')
hold on
plot(CG(5,:,1),Weight_CG(5,:),'k')
hold on
plot(CG(6,:,1),Weight_CG(6,:),'k')
hold on
plot(130.3132,Weightplot.We,'ro','linewidth',4)
hold on
plot(130.2768,Weightplot.Wo,'go','linewidth',4)
hold on
plot([aft,aft],[Weight_min-1e5,Weight_max+1e5], 'b--','linewidth',2)
hold on
plot([fwd,fwd],[Weight_min-1e5,Weight_max+1e5], 'b--','linewidth',2)
axis([fwd-0.1, aft+0.1, Weight_min-1e5, Weight_max+1e5])

set(gca,'FontSize',12)
xlabel('X_c_g Location [ft]')
ylabel('Weight [lbs]')
title('CG Excursion for Mk35: All Sequences')

H = legend('CGEmpty Weight', 'CG: Maximum Takeoff Weight','Location','NorthEast');
set(H,'FontSize',10)
end