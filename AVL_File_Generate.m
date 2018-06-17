function [] = AVL_File_Generate(Area, Main, Geom, tipTwist, rootTwist)
 
%AVL_File_Generate
%This file generates a .avl file corresponding to an airplane
%configuration. It is presented as we implemented in our code, for a
%lifting fuselage, single main wing, twin vertical tails, and H-tail with
%elevators for pitch moment balancing.

%The following code consists of two parts - a list of variables that you
%can plug in from your code; and a file generator corresponding to an AVL
%input file. It will output a new text input file in the current MATLAB
%directory

%% Reference values
Sref = Area.Wing; %wing area
Cref = Geom.Main.MAC; %reference chord??????????????????????????????????
Bref = Main.Span; %reference span
Xref = 35.794; %moment reference location????????????????????
Yref = -0.000; %moment reference location????????????????????
Zref = 3.776; %moment reference location???????????????????????
CDoref = 0.00; %This is Cd0 for cruise I MADE IT UP??????????????????????????????

%% Main Wing
%Root station
Xle_WingRoot = 0.0;
Yle_WingRoot = 0.0;
Zle_WingRoot = 0.0;
chord_WingRoot = Geom.Main.Chord_Root;
twist_WingRoot = 2;
Nspan_WingRoot = 19;
Sspace_WingRoot = 3;

sweep_main    = Geom.Main.sweep;
dihedral_main = Geom.Main.dihedral;

%Tip station
Yle_WingTip = Geom.Main.Span/2;
Xle_WingTip = Yle_WingTip*tand(sweep_main);
Zle_WingTip = Yle_WingTip*tand(dihedral_main);
chord_WingTip = Geom.Main.Chord_Root*Geom.Main.taper;
twist_WingTip = tipTwist;
Nspan_WingTip = 19;
Sspace_WingTip = 3;

sweep_tip    = Geom.Tip.sweep;
dihedral_tip = Geom.Tip.dihedral;

ANGLE_Wing = 0; %twist angle bias for whole surface

%% HARD-CODED VARIABLES (These are only defined here)
%% Main Wing
Nchord_Wing = 13;
Cspace_Wing = 1.0;
Nspan_Wing = 0; 
Sspace_Wing = 0;
YDUPLICATE_Wing = 0.0; %reflect image wing about y=0 plane
TRANSLATE_x_Wing = 0.0; %x bias for whole surface
TRANSLATE_y_Wing = 0.0; %y bias for whole surface
TRANSLATE_z_Wing = 0.0; %z bias for whole surface

%% OUTPUT - BEGIN WRITING TO .AVL FILE

% If the above variables don't make sense, replace the values below with your own inputs or hard-coded values. 
% It should correspond roughly to the sample .avl input file in the notes

%OPENS OUTPUT FILE
fileID = fopen('./AVL/airplane.avl','w');

%% WING
fprintf(fileID,'Spitfire_wing\r\n');
fprintf(fileID,'0.0                      | Mach\r\n');
fprintf(fileID,'0     0     0            | iYsym  iZsym  Zsym\r\n');
fprintf(fileID,'\r\n');

fprintf(fileID,'%7.5f %7.5f   %7.5f        Sref   Cref   Bref\r\n',...
    Sref,Cref,Bref);
fprintf(fileID,'# Note : check consistency of area unit above with length units of the file\r\n');
fprintf(fileID,'%7.5f %7.5f  %7.5f       Xref   Yref   Zref\r\n',...
    Xref,Yref,Zref);
fprintf(fileID,'%9.7f                    CDp  (optional)\r\n',CDoref);
fprintf(fileID,'\r\n');
fprintf(fileID,'\r\n');
fprintf(fileID,'\r\n');

fprintf(fileID,'#===================================================\r\n');
fprintf(fileID,'SURFACE             | (keyword)\r\n');
fprintf(fileID,'Test Plane_Wing\r\n');
fprintf(fileID,'#Nchord   Cspace  [ Nspan  Sspace ]\r\n');
fprintf(fileID,'%7.5f  %7.5f\r\n',Nchord_Wing,Cspace_Wing);
fprintf(fileID,'\r\n');

fprintf(fileID,'INDEX                | (keyword)\r\n');
fprintf(fileID,'22060                | Lsurf\r\n');
fprintf(fileID,'\r\n');

fprintf(fileID,'YDUPLICATE\r\n');
fprintf(fileID,'     %7.5f\r\n',YDUPLICATE_Wing);
fprintf(fileID,'\r\n');

fprintf(fileID,'SCALE\r\n');
fprintf(fileID,'1.0 1.0 1.0\r\n');
fprintf(fileID,'\r\n');

fprintf(fileID,'TRANSLATE\r\n');
fprintf(fileID,'    %7.5f     %7.5f     %7.5f\r\n',TRANSLATE_x_Wing,TRANSLATE_y_Wing,TRANSLATE_z_Wing);
fprintf(fileID,'\r\n');

fprintf(fileID,'ANGLE\r\n');
fprintf(fileID,'     %7.5f         | dAinc\r\n',ANGLE_Wing);
fprintf(fileID,'\r\n');
fprintf(fileID,'\r\n');

fprintf(fileID,'#First Section (Root)------------------------------------------------------\r\n');
fprintf(fileID,'#    Xle         Yle         Zle         chord       angle   Nspan  Sspace\r\n');
fprintf(fileID,'SECTION \r\n');
fprintf(fileID,'     %7.5f     %7.5f     %7.5f     %7.5f         %7.5f   %7.5f      %7.5f\r\n',...
    Xle_WingRoot,Yle_WingRoot,Zle_WingRoot,chord_WingRoot,twist_WingRoot,Nspan_WingRoot,Sspace_WingRoot);
fprintf(fileID,'AFIL\r\n');
fprintf(fileID,'.\\AVL\\SC(2)-0610.dat\r\n');
fprintf(fileID,'\r\n');
fprintf(fileID,'\r\n');

fprintf(fileID,'#Fourth Section (Tip)------------------------------------------------------\r\n');
fprintf(fileID,'SECTION \r\n');
fprintf(fileID,'     %7.5f      %7.5f     %7.5f     %7.5f       %7.5f   %7.5f     %7.5f\r\n',...
    Xle_WingTip,Yle_WingTip,Zle_WingTip,chord_WingTip,twist_WingTip,Nspan_WingTip,Sspace_WingTip);
fprintf(fileID,'AFIL\r\n');
fprintf(fileID,'.\\AVL\\SC(2)-0610.dat\r\n');
fprintf(fileID,'#--------------------------------------------------------------------------\r\n');

%CLOSES OUTPUT FILE
fclose(fileID);