function [] = AVL_CaseFile_Generate(Area, Main, W, CDo, CL)

%% Input values
Mach = Main.Mach;
vCru = Main.V_Cruise*1.688;
dens = Main.rho_c_rho_SL*1.225*0.00194;
mass = W*0.0310809502 ;

%% OUTPUT - BEGIN WRITING TO .AVL FILE
%OPENS OUTPUT FILE
fileID = fopen('./AVL/case','w');

%% WING
fprintf(fileID,'\r\n');
fprintf(fileID,'-----------------------------------------------------------\r\n');
fprintf(fileID,'Run case  1:  -SPITFIRE-\r\n');
fprintf(fileID,'\r\n');
fprintf(fileID,'alpha       -> CL    = %7.3f \r\n',CL);
fprintf(fileID,'beta        -> beta  = 0.00000\r\n');
fprintf(fileID,'pb/2V       -> pb/2V = 0.00000\r\n');
fprintf(fileID,'qc/2V       -> qc/2V = 0.00000\r\n');
fprintf(fileID,'rb/2V       -> rb/2V = 0.00000\r\n');
fprintf(fileID,'\r\n');
fprintf(fileID,'alpha     =  0.00000      deg\r\n');
fprintf(fileID,'beta      =  0.00000      deg\r\n');
fprintf(fileID,'pb/2V     =  0.00000\r\n');
fprintf(fileID,'qc/2V     =  0.00000\r\n');
fprintf(fileID,'rb/2V     =  0.00000\r\n');
fprintf(fileID,'CL        =  %7.3f \r\n',CL);
fprintf(fileID,'CDo       =  %7.3f \r\n',CDo);
fprintf(fileID,'bank      =  0.00000      deg\r\n');
fprintf(fileID,'elevation =  0.00000      deg\r\n');
fprintf(fileID,'heading   =  0.00000      deg\r\n');
fprintf(fileID,'Mach      =  %7.3f \r\n',Mach);
fprintf(fileID,'velocity  =  %7.3f      Lunit/Tunit\r\n',vCru);
fprintf(fileID,'density   =  %7.5f      Munit/Lunit^3\r\n',dens);
fprintf(fileID,'grav.acc. =  32.2000      Lunit/Tunit^2\r\n');
fprintf(fileID,'turn_rad. =  0.00000      Lunit \r\n');
fprintf(fileID,'load_fac. =  1.00000 \r\n');
fprintf(fileID,'X_cg      =  0.00000      Lunit\r\n');
fprintf(fileID,'Y_cg      =  0.00000      Lunit\r\n');
fprintf(fileID,'Z_cg      =  0.00000      Lunit\r\n');
fprintf(fileID,'mass      =  %7.1f      Munit\r\n',mass);
fprintf(fileID,'Ixx       =  1.00000      Munit-Lunit^2\r\n');
fprintf(fileID,'Iyy       =  1.00000      Munit-Lunit^2\r\n');
fprintf(fileID,'Izz       =  1.00000      Munit-Lunit^2\r\n');
fprintf(fileID,'Ixy       =  0.00000      Munit-Lunit^2\r\n');
fprintf(fileID,'Iyz       =  0.00000      Munit-Lunit^2\r\n');
fprintf(fileID,'Izx       =  0.00000      Munit-Lunit^2\r\n');
fprintf(fileID,'visc CL_a =  0.00000\r\n');
fprintf(fileID,'visc CL_u =  0.00000\r\n');
fprintf(fileID,'visc CM_a =  0.00000\r\n');
fprintf(fileID,'visc CM_u =  0.00000\r\n');

%CLOSES OUTPUT FILE
fclose(fileID);

end

