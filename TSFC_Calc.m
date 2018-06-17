function [TSFC] = TSFC_Calc(alt, M)

if M == 0
    TSFC = -6.4053e-7*alt + 0.29257;

elseif 0 < M && M <= 0.1
    c0 = -6.4053e-7*alt + 0.29257;
    c1 = -8.6960e-7*alt + 0.33025;
    r  = 0.1;    
    TSFC = (M-0)/r*(c1-c0)+c0;
  
elseif 0.1 < M && M <= 0.2
    c1 = -8.6960e-7*alt + 0.33025;
    c2 = -1.1163e-6*alt + 0.37009;
    r  = 0.1;    
    TSFC = (M-0.1)/r*(c2-c1)+c1;

elseif 0.2 < M && M <= 0.3
    c2 = -1.1163e-6*alt + 0.37009;
    c3 = -1.9619e-6*alt + 0.41442;
    r  = 0.1;    
    TSFC = (M-0.2)/r*(c3-c2)+c2;

elseif 0.3 < M && M <= 0.4
    c3 = -1.9619e-6*alt + 0.41442;
    c4 = -2.5310e-6*alt + 0.45862;
    r  = 0.1;    
    TSFC = (M-0.3)/r*(c4-c3)+c3;

elseif 0.4 < M && M <= 0.5
    c4 = -2.5310e-6*alt + 0.45862;
    c5 = -2.8540e-6*alt + 0.50875;
    r  = 0.1;    
    TSFC = (M-0.4)/r*(c5-c4)+c4;

elseif 0.5 < M && M <= 0.6
    c5 = -2.8540e-6*alt + 0.50875;
    c6 = 3.21e-15*alt^3 - 2.45e-10*alt^2 + 3.9e-6*alt + 0.488;
    r  = 0.1;    
    TSFC = (M-0.5)/r*(c6-c5)+c5;

elseif 0.6 < M && M <= 0.7
    c6 = 3.21e-15*alt^3 - 2.45e-10*alt^2 + 3.9e-6*alt + 0.488;
    c7 = 3.59e-15*alt^3 - 2.66e-10*alt^2 + 3.93e-6*alt + 0.534;
    r  = 0.1;    
    TSFC = (M-0.6)/r*(c7-c6)+c6;

elseif 0.7 < M && M <= 0.75
    c7 = 3.59e-15*alt^3 - 2.66e-10*alt^2 + 3.93e-6*alt + 0.534;
    c75 = 1.12e-10*alt^2 - 9.16e-6*alt + 0.702;
    r  = 0.05;    
    TSFC = (M-0.7)/r*(c75-c7)+c7;

elseif 0.75 < M && M <= 0.8
    c75 = 1.12e-10*alt^2 - 9.16e-6*alt + 0.702;
    c8 = 1.17e-10*alt^2 - 9.54e-6*alt + 0.728;
    r  = 0.05;    
    TSFC = (M-0.75)/r*(c8-c75)+c75;

elseif 0.8 < M && M <= 0.85
    c8 = 1.17e-10*alt^2 - 9.54e-6*alt + 0.728;
    c85 = 1.79e-10*alt^2 - 1.41e-5*alt + 0.832;
    r  = 0.05;    
    TSFC = (M-0.8)/r*(c85-c8)+c8;

elseif 0.85 < M && M <= 0.9
    c85 = 1.79e-10*alt^2 - 1.41e-5*alt + 0.832;
    c9 = 1e-10*alt^2 - 7.95e-6*alt + 0.731;
    r  = 0.05;    
    TSFC = (M-0.85)/r*(c9-c85)+c85;

else
    TSFC = -1.6250e-7*alt + 0.59197;

end

TSFC = TSFC*0.8;