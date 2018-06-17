function [P1,P2,P3] = T_S_Constraint(Req, Area, Main)

%% T_S Constraint Curves

s = linspace(2000,10e3,100);% Vector of Wing Area [ft^2]
t = linspace(0,5e5,100);% Vector of Thrust [lb]

% Determines Contraint curves for different flight parameters (takeoff, landing, cruise, etc.)
T_plot = zeros(size(s,2),10);
S_plot = zeros(size(t,2),1);

for j = 1:size(t,2)
    
    tolerance = 1e-3;
    
    for k = 1
        S = 10000;
        dif = 1;
        while dif > tolerance

            S_old           = S;
            Area.Wing       = S_old;
            Area.Wing_CS    = Area.Wing*0.3;
            Area.Wing_no_CS = Area.Wing*0.7;
            Main.Wing_Area  = S_old;
            Main.Thrust     = t(1,j)/2;
            Main.AR         = (Main.Span^2)/Main.Wing_Area;% Aspect Ratio      
            [Weight]        = Spitfire_Weight(Req, Area, Main);
            W               = Weight.MTOW;
            
            W_S             = W/S_old;
            [TW, WS]        = Spitfire_TS_TW(Main, Weight, k, W_S);
            S               = W/WS;
            
            dif             = abs(S - S_old);

        end
        S_plot(j,k) = S;
    end

end

for j = 1:size(s,2)
    
    tolerance = 1e-3;
    
    for k = 2:7
        T = 1e5;
        dif = 1;
        while dif > tolerance

            T_old           = T;
            Main.Thrust     = T_old/2;
            Main.Wing_Area  = s(1,j);
            Area.Wing       = s(1,j);
            Area.Wing_CS    = Area.Wing*0.3;
            Area.Wing_no_CS = Area.Wing*0.7;
            Main.AR               = (Main.Span^2)/Main.Wing_Area;% Aspect Ratio
            [Weight]        = Spitfire_Weight(Req, Area, Main);
            W               = Weight.MTOW;
            W_S             = W/s(1,j);
            [TW, WS]        = Spitfire_TS_TW(Main, Weight, k, W_S);
            T               = TW*W;
            dif             = abs(T - T_old);

        end
        T_plot(j,k) = T;
    end
    
end

S_min = zeros(size(s,2),1);
rho   = Main.rho_rho_SL*1.225*0.00194;%[slug/ft^3]
for j = 1:size(t,2)
    
    for k = 1:size(s,2)
        
        Main.Thrust     = t(1,j)/2;
        Main.Wing_Area  = s(1,k);
        Area.Wing       = s(1,k);
        Area.Wing_CS    = Area.Wing*0.3;
        Area.Wing_no_CS = Area.Wing*0.7;
        Main.AR               = (Main.Span^2)/Main.Wing_Area;% Aspect Ratio
        [Weight]        = Spitfire_Weight(Req, Area, Main);
        V_stall         = sqrt(2*Weight.MTOW/(rho*Main.Wing_Area*Main.CL_max));
        S_min(j,k)      = 2*Weight.MTOW/(rho*0.9*(V_stall + 20*1.68781)^2);
      
    end
    
end

P1 = zeros(1,2);
P2 = zeros(6,3);
P3 = zeros(1,2);

P1(1,:) = polyfit(t,S_plot(:,1)',1);% s = A*t + B
S1      = polyval(P1(1,:),t);

P2(1,:) = polyfit(s,T_plot(:,2)',2);% t = A*s^2 + B*s + c
S2      = polyval(P2(1,:),s);

P2(2,:) = polyfit(s,T_plot(:,3)',2);% t = A*s^2 + B*s + c
S3      = polyval(P2(2,:),s);

P2(3,:) = polyfit(s,T_plot(:,4)',2);% t = A*s^2 + B*s + c
S4      = polyval(P2(3,:),s);

P2(4,:) = polyfit(s,T_plot(:,5)',2);% t = A*s^2 + B*s + c
S5      = polyval(P2(4,:),s);

P2(5,:) = polyfit(s,T_plot(:,6)',2);% t = A*s^2 + B*s + c
S6      = polyval(P2(5,:),s);

P2(6,:) = polyfit(s,T_plot(:,7)',2);% t = A*s^2 + B*s + c
S7      = polyval(P2(6,:),s);

P3(1,:) = polyfit(t,S_min(:,1)',1);% s = A*t^2 + B*t + c
S8      = polyval(P3(1,:),t);

figure(1)
clf
plot(S1,t, 'LineWidth', 2, 'color', [1 0.4 0])
hold on
plot(s,S2, 'LineWidth', 2, 'color', [1 0 0])
hold on
plot(s,S3, 'LineWidth', 2, 'color', [1 1 0])
hold on
plot(s,S4, 'LineWidth', 2, 'color', [0 1 1])
hold on
plot(s,S5, 'LineWidth', 2, 'color', [0 1 0])
hold on
plot(s,S6, 'LineWidth', 2, 'color', [1 0 1])
hold on
plot(s,S7, 'LineWidth', 2, 'color', [0.8 0.8 0.8])
hold on
plot(S8,t,'--', 'LineWidth', 2, 'color','k')
hold on

legend('Landing','Takeoff','1st Climb','2nd Climb','3rd Climb',...
       'Ceiling','Cruise','Location','EastOutside')
xlabel('Wing Area (S) [ft^2]','fontsize',16)
ylabel('Thrust (T) [lbf]','fontsize',16)
set(gca,'FontSize',16)
axis([2000, 10000, 0, 5e5])
end
