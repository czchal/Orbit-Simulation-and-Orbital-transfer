% Define Parameters
param.mu = 398601.2;          %[km^3/s^2]
param.Re = 6378;              %[km] 
param.cd= 2.07;             % cd for ISS
param.rho= 3.8*10^-12;       %  [kg m-3] mean density at the altitude of ISS 
param.S= 1.5*10^-3;           %since the surface area of ISS changes from 700m^2 to 23000m^2 an average of the two values.  
param.m = 250000;                 %[kg]   mass of ISS
pi=3.14159265359;

% for inital condition of r and v calculate the required parameters like a, T and dt
% 
r_eci= [3584.681724590060 -2531.834101518260 5181.093477824530]; 
v_eci= [5.55481871603161 5.10031844394557 -1.35070548166970];



%% Simulate orbit for 5 period solving two-body dynamics using RK-4

X=num_orbit_nodrag(r_eci,v_eci,5,param);

figure(2);
plot3(X(:,1),X(:,2),X(:,3),'r',LineWidth=2)   %% plot for Problem I-6
xlabel('ECI x [km]');
ylabel('ECI y [km]');
zlabel('ECI z [km]');
title('Satellite Orbit in ECI Coordinates');
grid on
%% Simulate orbit for 5 period solving two-body dynamics using RK-4 considering drag due to atmosphere

[X,semi_major]=num_orbit_withdrag(r_eci,v_eci,50,param);

figure(3);
plot3(X(:,1),X(:,2),X(:,3),'r',LineWidth=2)   %% plot for Problem I-7
xlabel('ECI x [km]');
ylabel('ECI y [km]');
zlabel('ECI z [km]');
title('Satellite Orbit in ECI Coordinates');
grid on
figure(4);
plot(1:1:length(semi_major),semi_major)
