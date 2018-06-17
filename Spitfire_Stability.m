function [ SM, xnp, SM_empty ] = Spitfire_Stability(X,x_cg, x_cg_empty)
[Req, Area, Main, Geom] = Variables(X);
Weight = Spitfire_Weight(Req, Area, Main, Geom, X);

%% Inputs

WTO   = [Weight.We Weight.MTOW];                                           % Empty weight and takeoff weight (lbs)
taper = Geom.Main.taper;                                                   % Taper ratio
b     = Geom.Main.Span;                                                    % Wingspan
x_hs  = Geom.HT.Dist;                                                      % Distance from nose to front of horizontal stabilizer (ft)

%% Wing

y_bar      = Geom.Main.Y_bar;                                              % Spanwise location of MAC (ft)
sweep      = Geom.Main.sweep;                                              % Sweep angle in degrees
AR         = Geom.Main.AR;                                                 % Aspect ratio
Sw         = Area.Wing;                                                    % Wing reference area (ft^2)
c          = Geom.Main.MAC;                                                % Mean aerodynamic chord (ft)
eta        = 0.97;                                                         % Section Lift Curve Difference Coefficient
CLaw       = 2*pi*AR/(2+sqrt((AR/eta)^2*(1+(tand(sweep))^2 - X(3)^2) + 4));  
xmac_wing  = Geom.Main.Wing_Dist + tand(sweep)*y_bar + c/4;                % Distance from nose to 1/4 MAC (ft)

%% Horizontal stabilizer

Sh       = Area.HT;                                                        % Area of hoizontal stabilizer (ft^2)
AR_HT    = Geom.HT.AR;                                                     % Aspect ratio of Horizontal Stabilizer
sweep_HT = Geom.HT.sweep;
CLah0    = 2*pi*AR_HT/(2+sqrt((AR_HT/eta)^2*...
          (1+(tand(sweep_HT))^2 - X(3)^2) + 4));
xmac_HT  = Geom.HT.Dist + tand(sweep_HT)*Geom.HT.Y_bar + Geom.HT.MAC/4;    % Distance from nose to 1/4 MAC (ft)
lh       = xmac_HT - xmac_wing;                                            % Distance from wing front edge to horizontal stabilizer (ft)
de_da    = 2*CLaw/(pi*AR);                                                 % Downwash effect slope
CLah     = CLah0*(1-de_da)*0.9;                                            % CL vs alpha slope for horizontal stabilizer, accounting for main wing downwash

%% Fuselage

lf = Geom.Fuse.length;                                                     % Length of fuselage                                             
wf = Geom.Fuse.Diam;                                                       % Maximum width of fuselage (ft)
X  = xmac_wing/lf;                                                         % Wing 1/4 chord position reletaive to aircraft length
Kf = 12.0076*X^4 - 20.0177*X^3 + 12.8458*X^2 - 1.9468*X + 0.1980;          % Empirical factor
dCmfus_dCl = Kf*wf.^2*lf/(Sw*c*CLaw);                                      % Moment derivative for fuselage

%% Stability analysis

xcg = x_cg - xmac_wing;                                                    % Position of CG relative to 1/4 MAC (ft)
xnp = c * (lh*Sh*CLah/(c*Sw*CLaw) - dCmfus_dCl);
SM = (xnp-xcg)/c;
xcg_empty = x_cg_empty - xmac_wing;                                        % Position of CG relative to 1/4 MAC (ft)
SM_empty = (xnp-xcg_empty)/c;

end