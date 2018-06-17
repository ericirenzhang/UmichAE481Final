function Obj = ObjectiveFunction(X)
Scale = [5000 200000 0.84 35000 30 0.2 100];
X = X.*Scale;
disp(['Current Values: ',num2str(X(1)), ', ',num2str(X(2)),...
      ', ',num2str(X(3)),', ',num2str(X(4)),', ',num2str(X(5)),...
      ', ',num2str(X(6)),' and ',num2str(X(7))]);

%% Objective Function
[Req, Area, Main, Geom] = Variables(X);

PlanePlot(X);

% Other Variables
Ranges    = [8166,7366,4427,3691,1806];
Weighting = [1 2 6 8 1];
ObjFun    = 0;

for i = 1:5
            
    Req.Airport_Range = Ranges(i);

    [Weight, Flight] = Spitfire_Weight(Req, Area, Main, Geom, X);
    [L_D] = lift_over_drag(Area, Main, Geom, Weight);
    Main.L_D = L_D;
    Main.Airport_t_b = ceil((Req.Airport_Range/Main.V_Cruise + 1)*2)/2;

    [Cost] = CostCOC(Weight, Req, Main);
    ObjFun = ObjFun + Cost.COC/Req.Airport_Range*Weighting(i);

end

Obj = 1/(18*Req.Passengers)*ObjFun*100;

end