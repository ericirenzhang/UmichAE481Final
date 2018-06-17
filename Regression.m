%%  (Thrust, SFC) vs. (Mach, Altitude)

clear all; close all; clc;

A = textread('GE-90.txt');

Mach = A(:,1);
Altitude = A(:,2);
powerCode = A(:,3);
T_gross = A(:,4);
fuelBurn = A(:,6);
TSFC = A(:,7);

C = unique(powerCode);

M0  = [];
M1  = [];
M2  = [];
M3  = [];
M4  = [];
M5  = [];
M6  = [];
M7  = [];
M75 = [];
M8  = [];
M85 = [];
M9  = [];



for i = 1:length(Mach)

    if Mach(i) < 0.4
        if powerCode(i) == C(7)
            if Mach(i) == 0
                M0 = [M0;Altitude(i), TSFC(i), T_gross(i)];
            elseif Mach(i) == 0.1
                M1 = [M1;Altitude(i), TSFC(i), T_gross(i)];
            elseif Mach(i) == 0.2
                M2 = [M2;Altitude(i), TSFC(i), T_gross(i)];
            else
                M3 = [M3;Altitude(i), TSFC(i), T_gross(i)];
            end
        end
    
    
    elseif Mach(i) == 0.4
        if powerCode(i) == C(6)
            M4 = [M4;Altitude(i), TSFC(i), T_gross(i)];
        end
        
    elseif Mach(i) == 0.5
        if powerCode(i) == C(5)
            M5 = [M5;Altitude(i), TSFC(i), T_gross(i)];
        end
        
    else
        if powerCode(i) == C(7)
            if Mach(i) == 0.6
                M6 = [M6;Altitude(i), TSFC(i), T_gross(i)];
            end
            if Mach(i) == 0.7
                M7 = [M7;Altitude(i), TSFC(i), T_gross(i)];
            end
            if Mach(i) == 0.75
                M75 = [M75;Altitude(i), TSFC(i), T_gross(i)];
            end
            if Mach(i) == 0.8
                M8 = [M8;Altitude(i), TSFC(i), T_gross(i)];
            end
            if Mach(i) == 0.85
                M85 = [M85;Altitude(i), TSFC(i), T_gross(i)];
            end
            if Mach(i) == 0.9
                M9 = [M9;Altitude(i), TSFC(i), T_gross(i)];
            end
        end        
    end
end

figure()
plot(M0(:,1),M0(:,2),'.',M1(:,1),M1(:,2),'.',...
    M2(:,1),M2(:,2),'.',M3(:,1),M3(:,2),'.')
xlabel('Altitude [ft]')
ylabel('TSFC [lb/(lbf-hr)????]')
title('100% Throttle, 0 < M < 0.4')
axis([0, 1000 + max([M0(:,1);M1(:,1);M2(:,1);M3(:,1)]),...
    0, .05 + max([M0(:,2);M1(:,2);M2(:,2);M3(:,2)])])
legend('M = 0.0','M = 0.1','M = 0.2','M = 0.3')

figure()
plot(M4(:,1),M4(:,2),'.')
xlabel('Altitude [ft]')
ylabel('TSFC [lb/(lbf-hr)????]')
title('90% Throttle, M = 0.4')
axis([0, 1000 + max(M4(:,1)), 0, .05 + max(M4(:,2))])
legend('M = 0.4')

figure()
plot(M5(:,1),M5(:,2),'.')
xlabel('Altitude [ft]')
ylabel('TSFC [lb/(lbf-hr)????]')
title('80% Throttle, M = 0.5')
axis([0, 1000 + max(M5(:,1)), 0, .05 + max(M5(:,2))])
legend('M = 0.5')

figure()
plot(M6(:,1),M6(:,2),'.',M7(:,1),M7(:,2),'.',...
    M75(:,1),M75(:,2),'.',M8(:,1),M8(:,2),'.',...
    M85(:,1),M85(:,2),'.',M9(:,1),M9(:,2),'.')
xlabel('Altitude [ft]')
ylabel('TSFC [lb/(lbf-hr)????]')
title('100% Throttle, 0.6 < M < 0.9')
axis([0, 1000 + max([M6(:,1);M7(:,1);M75(:,1);M8(:,1);M85(:,1);M9(:,1)]),...
    0, .05 + max([M6(:,2);M7(:,2);M75(:,2);M8(:,2);M85(:,2);M9(:,2)])])
legend('M = 0.6','M = 0.7','M = 0.75','M = 0.8','M = 0.85','M = 0.9')
      
            
            
            
            