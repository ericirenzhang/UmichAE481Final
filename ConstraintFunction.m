function [c,ceq] = ConstraintFunction(X)
Scale = [5000 200000 0.84 35000 30 0.2 100];
X = X.*Scale;
%% Nonlinear inequality constraints
[Req, Area, Main, Geom] = Variables(X);
[Weight, Flight] = Spitfire_Weight(Req, Area, Main, Geom, X);
[L_D, C_L, C_D] = lift_over_drag(Area, Main, Geom, Weight);
Main.L_D = L_D;
[C, C2] = Spitfire_TS_TW(Main,C_D, C_L);
W = Weight.MTOW;
[X_CG_TO, X_CG_Empty, x_LG] = Center_of_gravity_new(Geom, Weight);
[ SM, X_NP, SM_empty ] = Spitfire_Stability(X,X_CG_TO,X_CG_Empty);

c = [W/X(1)*C(1,1) + X(2)/W*C(1,2) + C(1,3);
     W/X(1)*C(2,1) + X(2)/W*C(2,2) + C(2,3);
     W/X(1)*C(3,1) + X(2)/W*C(3,2) + C(3,3);
     W/X(1)*C(4,1) + X(2)/W*C(4,2) + C(4,3);
     W/X(1)*C(5,1) + X(2)/W*C(5,2) + C(5,3);
     W/X(1)*C(6,1) + X(2)/W*C(6,2) + C(6,3);
     X(1)/W*C2(1,1) + W/X(1)*C2(1,2) + X(2)/W*C2(1,3);
     8  - atand((Geom.LG.Z_pos_fuse)/((2039 + 400)/12 - Geom.LG.X_pos))
     0.05 - SM;
     0.05 - SM_empty;
     SM - 0.40;
     SM_empty - 0.40;
     10 - atand((Geom.LG.X_pos - X_CG_TO)/Geom.LG.Z_pos);
     10 - atand((Geom.LG.X_pos - X_CG_Empty)/Geom.LG.Z_pos);
     40 - Geom.Main.Chord_Root;
     Flight.Ceil - 47000];

%% Nonlinear equality constraints

ceq = [];