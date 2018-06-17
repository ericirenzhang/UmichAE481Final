function [ObjFun] = T_S_Plot(Req, Area, Main, Geom, X)
Ranges = [8166,7366,4427,3691,1806];
Weighting = [1 2 6 8 1];

s = linspace(4000,10e3,100);% Vector of Wing Area [ft^2]
t = linspace(30000,5e5,100);% Vector of Thrust [ft^2]

COC = zeros(length(s),length(t)); %empty matrix
ObjFun = zeros(length(s),length(t)); %empty matrix
L_D_temp = zeros(length(s),length(t)); %empty matrix

% Determines COC with varying thrust and wing area
for i = 1:5
    for j = 1:size(s,2)

        for k = 1:size(t,2)

            [Req, Area, Main, Geom] = Variables([s(1,j) t(1,k) X(3:end)]); 
            Req.Airport_Range = Ranges(i);
            [Weight, Flight] = Spitfire_Weight(Req, Area, Main, Geom, X);
            [L_D, C_L, C_D] = lift_over_drag(Area, Main, Geom, Weight);
            Main.L_D                = L_D;
            L_D_temp(j,k)           = Main.L_D;
            Main.Airport_t_b = ceil((Req.Airport_Range/Main.V_Cruise + 1)*2)/2;
            [Cost]                  = CostCOC(Weight, Req, Main);
            ObjFun(j,k)             = ObjFun(j,k) + Cost.COC/Req.Airport_Range*Weighting(i);
            
        end

    end
end
ObjFun = 1/(18*Req.Passengers).*ObjFun*100;

% Plots Contour of COC
figure(2)
[S,T] = meshgrid(s,t);
contourf(S,T,ObjFun',25); 
c = colorbar('v');
set(get(c,'ylabel'),'string','Objective Function [Cents/pax-mile]');

%% T_S Constraint Curves

s_full = linspace(4000,10e3,100);% Vector of Wing Area [ft^2]
t_full = linspace(30000,5e5,100);% Vector of Thrust [lb]
s = s_full'*ones(1,6);
t = t_full'*ones(1,1);

% Determines Contraint curves for different flight parameters (takeoff, landing, cruise, etc.)
T_plot = zeros(size(s_full,2),6);
S_plot = zeros(size(t_full,2),1);

for j = 1:size(t,1)
    
    tol = 1e-3;
    
    for k = 1
        S = 10000;
        i = 1;
        conv = 0;
        
        while conv == 0

            if i > 200
                S_plot(j,k) = 0;
                t(j,k) = 0;
                warning('Failed to Converge');
                break
            end
            
            S_old                   = S;
            [Req, Area, Main, Geom] = Variables([S_old t(j,1) X(3:end)]);     
            [Weight, Flight] = Spitfire_Weight(Req, Area, Main, Geom, X);
            [L_D, C_L, C_D] = lift_over_drag(Area, Main, Geom, Weight);
            Main.L_D                = L_D;         
            C_L_max                 = C_L + 1.56;            
            S = Weight.MTOW/(((Main.rho_rho_SL)/80)*(Main.S_Land - 1000)/1.67*C_L_max/0.75);
            
            if abs(S - S_old) <= tol
                S_plot(j,k) = S;
                conv = 1;
            end

            i = i + 1;

        end
    end

end

for j = 1:size(s,1)
    
    tol = 1e-3;
    
    for k = 2:7
        T = 1.5e5;
        i = 1;
        conv = 0;
        while conv == 0

            if i > 200
                T_plot(j,k-1) = 0;
                s(j,k-1) = 0;
                warning('Failed to Converge');
                break
            end
            
            T_old                   = T;
            [Req, Area, Main, Geom] = Variables([s(j,1) T_old X(3:end)]); 
            [Weight, Flight] = Spitfire_Weight(Req, Area, Main, Geom, X);
            [L_D, C_L, C_D] = lift_over_drag(Area, Main, Geom, Weight);
            Main.L_D                = L_D;
            W                       = Weight.MTOW;
            W_S                     = W/s(j,1);
            [TW, WS]                = T_S_Plot_constraint(Main, Weight,k,W_S, Area, C_L, C_D);
            T                       = TW*W;
            
            if abs(T - T_old) <= tol
                T_plot(j,k-1) = T;
                conv = 1;
            end            
            
            i = i + 1;

        end
    end
    
end

S_plot_new = sparse(S_plot);
T_plot_new = sparse(T_plot);
s_new = sparse(s);
t_new = sparse(t);

S_min = zeros(size(s_full,2),1);
rho   = Main.rho_rho_SL*1.225*0.00194;%[slug/ft^3]

figure(2)
hold on
plot(S_plot_new(:,1),t_new(:,1), 'LineWidth', 2, 'color', [1 0.4 0]);% Landing
hold on;
plot(s_new(:,1),T_plot_new(:,1), 'LineWidth', 2, 'color', [1 0 0]);% Takeoff
hold on;
plot(s_new(:,2),T_plot_new(:,2), 'LineWidth', 2, 'color', [1 1 0]);%1st Climb
hold on;
plot(s_new(:,3),T_plot_new(:,3), 'LineWidth', 2, 'color', [0 1 1]);%2nd Climb
hold on;
plot(s_new(:,4),T_plot_new(:,4), 'LineWidth', 2, 'color', [0 1 0]);%3rd Climb
hold on;
plot(s_new(:,5),T_plot_new(:,5), 'LineWidth', 2, 'color', [1 0 1]);%Ceiling
hold on;
plot(s_new(:,6),T_plot_new(:,6), 'LineWidth', 2, 'color', [0.8 0.8 0.8]);%Cruise
hold on
[Req, Area, Main, Geom] = Variables(X);
Area_min = Geom.Main.Span*40*(1+Geom.Main.taper)/2;
plot([Area_min,Area_min],[30000,5e5],'--', 'LineWidth', 2, 'color','k')

legend('Objective Function','Landing','Takeoff','1st Climb','2nd Climb','3rd Climb',...
       'Ceiling','Cruise','Geom Constraint','Location','EastOutside')
xlabel('Wing Area (S) [ft^2]','fontsize',16)
ylabel('Thrust (T) [lbf]','fontsize',16)
set(gca,'FontSize',16)
axis([4000, 10000, 30000, 5e5])

end
