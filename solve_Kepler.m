function [E] = solve_Kepler(e,M)
%SOLVE_KEPLER solves kepler equation to obtain the value of Eccentric
%anolmay given eccentricy and true anomaly.Kepler's equation => M= E-esinE
E=M;
delta_M=E-e*sin(E)-M;
delta_E=1;
while abs(delta_E)>1e-4
    delta_E=-delta_M/(1-e*cos(E));
    E=E+delta_E;
    delta_M=E-e*sin(E)-M;
end
end

