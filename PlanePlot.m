function [] = PlanePlot(X)
%% Variables
[Req, Area, Main, Geom] = Variables(X);
[Weight, Flight] = Spitfire_Weight(Req, Area, Main, Geom, X);
[X_CG_TO, X_CG_Empty, x_LG] = Center_of_gravity_new(Geom, Weight);
[ SM, X_NP, SM_empty ] = Spitfire_Stability(X,X_CG_TO,X_CG_Empty);

%% Plot Wing
C_Root = Geom.Main.Chord_Root;

Ax = 0;
Bx = Geom.Main.Span/2*tand(Geom.Main.sweep);
Cx = Geom.Main.Span/2*tand(Geom.Main.sweep) + C_Root*Geom.Main.taper;
Dx = C_Root;
Ex = Geom.Main.Span/2*tand(Geom.Main.sweep) + C_Root*Geom.Main.taper;
Fx = Geom.Main.Span/2*tand(Geom.Main.sweep);
Gx = 0;

Ay = 0;
By = -Geom.Main.Span/2;
Cy = -Geom.Main.Span/2;
Dy = 0;
Ey = Geom.Main.Span/2;
Fy = Geom.Main.Span/2;
Gy = 0;

MAC       = (Geom.Main.MAC);%Wing Mean Aerodynamic Chord
Y_bar     = (Geom.Main.Y_bar);%Wing Mean Aerodynamic Chord Y location
MAC_LE    = Y_bar*tand(Geom.Main.sweep);%Leading Edge Position of Mean Aerodynamic Chord
Wing_Dist = Geom.Main.Wing_Dist;
MAC_quat  = Wing_Dist + Y_bar*tand(Geom.Main.sweep) + 0.25*MAC;

%% Plot Horizontal Tail
C_Root = Geom.HT.Chord_Root;

A1x = 0;
B1x = Geom.HT.Span/2*tand(Geom.HT.sweep);
C1x = Geom.HT.Span/2*tand(Geom.HT.sweep) + C_Root*Geom.HT.Taper;
D1x = C_Root;
E1x = Geom.HT.Span/2*tand(Geom.HT.sweep) + C_Root*Geom.HT.Taper;
F1x = Geom.HT.Span/2*tand(Geom.HT.sweep);
G1x = 0;

A1y = 0;
B1y = -Geom.HT.Span/2;
C1y = -Geom.HT.Span/2;
D1y = 0;
E1y = Geom.HT.Span/2;
F1y = Geom.HT.Span/2;
G1y = 0;

%% Plot Fuselage
nose_minor = 125;
nose_major = 400;
nose_x     = linspace(0,nose_major,100);
nose_y     = sqrt((nose_minor^2)*(1-(nose_x.^2)/nose_major^2));

tail_minor = 125;
tail_major = 750;
tail_x     = linspace(0,tail_major,100);
tail_y     = sqrt((tail_minor^2)*(1-(tail_x.^2)/tail_major^2));

%% Plot
figure(10)
clf
set(0,'DefaultLineLinewidth',2);
plot([MAC_LE (MAC_LE + MAC)] + Wing_Dist,[Y_bar Y_bar],'r-')
hold on
plot(X_CG_TO,0,'ro','linewidth',5)
hold on
plot(X_CG_Empty,0,'go','linewidth',5)
hold on
plot(MAC_quat+X_NP,0,'ko','linewidth',5)
hold on
plot([MAC_LE (MAC_LE + MAC)] + Wing_Dist,[-Y_bar -Y_bar],'r-')
hold on
plot([Ax Bx Cx Dx Ex Fx Gx] + Wing_Dist,[Ay By Cy Dy Ey Fy Gy],'b-') 
hold on
plot([A1x B1x C1x D1x E1x F1x G1x] + Geom.Fuse.length*0.86,[A1y B1y C1y D1y E1y F1y G1y],'b-') 
hold on
plot((-nose_x+400)./12,nose_y./12,'b-')
hold on
plot((-nose_x+400)./12, -nose_y./12,'b-')
hold on
plot([400 2439]./12,[-125 -125]./12,'b-')
hold on
plot((2439 + tail_x)./12,-tail_y./12,'b-')
hold on
plot((2439 + tail_x)./12,tail_y./12,'b-')
hold on
plot([2439 400]./12,[125 125]./12,'b-')

axis square
axis([-100 300, -200 200])
ylabel('Span [ft]','FontSize',16)
xlabel('Chord [ft]','FontSize',16)
title({'Wing Optimization',['SM (MTOW): ',num2str(SM),'    SM (Wempty): ',num2str(SM_empty)],['MTOW: ',num2str(Weight.MTOW),'    Wing Area: ',num2str(X(1))]},...
      'FontSize',16)
legend(['MAC: ', num2str(MAC),' [ft]'],['CG (MTOW): ',num2str(X_CG_TO),' [ft]'],...
       ['CG (Wempty): ',num2str(X_CG_Empty),' [ft]'],['NP: ',num2str(MAC_quat+X_NP),' [ft]'],'location','NorthWest')

end

