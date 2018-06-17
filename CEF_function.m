function [ CEF ] = CEF_function( b_year, t_year )
%This Function Calculates the Cost Escalation Factor (CEF)

% Equations
b_CEF = 5.17053 + 0.104981*(b_year - 2006);% Base CEF
t_CEF = 5.17053 + 0.104981*(t_year - 2006);% Then CEF

CEF = t_CEF/b_CEF;% Cost Escalation Factor (CEF)

end