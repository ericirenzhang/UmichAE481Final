function [AVL] = AVLrun(X,W,CDo,CL)

[~, Area, Main, Geom] = Variables(X);
AVL_CaseFile_Generate(Area, Main, W,CDo,CL);

tipAngle = 0;
rootAngle = 0;

AVL = zeros(4,1);

AVL_File_Generate(Area, Main, Geom, tipAngle, rootAngle);
[AVLdata] = runAVL();

AVL(1) = AVLdata.Alpha;
AVL(2) = AVLdata.CLtot;
AVL(3) = AVLdata.CDtot;
AVL(4) = AVLdata.e;

        