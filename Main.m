%% This script is written to simulate orbit solving Kepler's equation 

clc;clear;
%% Problem I
%initial values in KM and KM/s
r_eci= [3584.681724590060 -2531.834101518260 5181.093477824530]; 
v_eci= [5.55481871603161 5.10031844394557 -1.35070548166970];
% Define Parameters
param.mu = 398601.2;          %[km^3/s^2]
param.Re = 6378;              %[km] 
param.cd=  1.80;             % cd for ISS
param.rho= 3.8*10^-3;        % [kg km-3] mean density at the altitude of ISS 
param.S=   1.96462*10^-3;           %since the surface area of ISS changes from 700m^2 to 2300m^2 an average of the two values.  
param.m =  450000;              %[kg]   mass of ISS
pi=3.14159265359;
%% Problem I-1
[a,e,i,omega,argument_of_perigee,true_anomaly]=ECI2classical(r_eci,v_eci)

%% Problem I-2
T=2*pi()*a^1.5/param.mu^0.5
%% Problem I-4: track true anomaly change for 5 orbital period solving Kepler's equation using Newton Rapson method
t_a=track_true_anomaly(r_eci,v_eci,5);

figure(1)
plot(1:length(t_a),t_a)    % plot for problem I-4
xlabel('ECI t [seconds]');
ylabel('ECI true anomaly [rads]');
title('True anomaly trend over time');
%% Problem I-5:  simulate the trajectory of satellite in cartesian coordinate 
[positions,velocities]=simulate_orbit(r_eci,v_eci,5);

% plot for Problem I-5
figure(2);
plot3(positions(:,1),positions(:,2),positions(:,3),'r',LineWidth=2)   
xlabel('ECI x [km]');
ylabel('ECI y [km]');
zlabel('ECI z [km]');
title('Satellite Orbit in ECI Coordinates solving Keplers equation');
grid on

%% Problem I-6: simulate the trajectory using numerical intergration 

X=num_orbit_nodrag(r_eci,v_eci,5,param);

figure(3);
plot3(X(:,1),X(:,2),X(:,3),'r',LineWidth=2)   %% plot for Problem I-6
xlabel('ECI x [km]');
ylabel('ECI y [km]');
zlabel('ECI z [km]');
title('Satellite Orbit in ECI Coordinates using numerical integration');
grid on
%% Problem I-7: Simulate orbit for number of periods solving two-body dynamics using RK-4 considering drag due to atmosphere

% r_eci=[6021.49496554112	-4252.26846388357	8702.09881876200];
% v_eci=[2.76554969969890	1.88379621307132	0.00899440252689820];
% 
r_eci= [3584.681724590060 -2531.834101518260 5181.093477824530]; 
v_eci= [5.55481871603161 5.10031844394557 -1.35070548166970];

[X,semi_major]=num_orbit_withdrag(r_eci,v_eci,500,param);

figure(4);
plot3(X(:,1),X(:,2),X(:,3),'r',LineWidth=2)   %% plot for Problem I-7
xlabel('ECI x [km]');
ylabel('ECI y [km]');
zlabel('ECI z [km]');
title('Satellite Orbit in ECI Coordinates numerical integration with drag');
grid on
figure(5);
plot(1:1:length(semi_major),semi_major)
xlabel('time [seconds]');
ylabel('Semi major axis [km]');
title('Semi major axis change due to drag');
%% Probelm II 

%initial values in KM and KM/s

r_eci= [6000 5000 3000]; 
v_eci= [6.802 6.802 0];

mu= 398601.2; 
pi=3.14159265359;

% Define Parameters
param.mu = 398601.2;          %[km^3/s^2]
param.Re = 6378;              %[km] 
% param.cd=  2.07;             % cd for ISS
% param.rho= 3.8*10^-3;       %  [kg m-3] mean density at the altitude of ISS 
% param.S=  1.5*10^-3;           %since the surface area of ISS changes from 700m^2 to 23000m^2 an average of the two values.  
% param.m = 450000;                 %[kg]   mass of ISS
pi=3.14159265359;

%% Problem II-1-1

H=cross(r_eci,v_eci);
p=norm(H)^2/mu;

[a,e,i,omega,argument_of_perigee,true_anomaly]=ECI2classical(r_eci,v_eci)    %ECI to classical orbital elements convertion Problem I-1

%% Problem II-1-2
T=2*pi()*a^1.5/mu^0.5

%% Problem II-1-4
t_a=track_true_anomaly(r_eci,v_eci,5);

% plot for problem I-4
figure(6)
plot(1:length(t_a),t_a)  
xlabel('ECI t [seconds]');
ylabel('ECI true anomaly [rads]');

title('True anomaly trend over time');

%% Problem II-1-5
[positions,velocities]=simulate_orbit(r_eci,v_eci,1);

% plot for Problem I-5
figure(7);
plot3(positions(:,1),positions(:,2),positions(:,3),'r',LineWidth=2)   
xlabel('ECI x [km]');
ylabel('ECI y [km]');
zlabel('ECI z [km]');
title('Satellite Orbit in ECI Coordinates');
grid on

%% Problem II-1-6: simulate the trajectory using numerical intergration 

X=num_orbit_nodrag(r_eci,v_eci,5,param);

figure(8);
plot3(X(:,1),X(:,2),X(:,3),'r',LineWidth=2)   %% plot for Problem I-6
xlabel('ECI x [km]');
ylabel('ECI y [km]');
zlabel('ECI z [km]');
title('Satellite Orbit in ECI Coordinates using numerical integration');
grid on
%% Problem III: Simulate orbit transfer between two circular orbit
param.mu = 398601.2;          %[km^3/s^2]
param.Re = 6378;  
pi=3.14159265359;
mu=398601.2;
a1=param.Re+200;
% a2=42167;
a2=a1*15;


% generate the initial position and velocity 
e1=0;
i1=0;
omega1=0;
argument_of_perigee1=90*pi/180;
true_anomaly1=0;

% [a,e,i,omega,argument_of_perigee,true_anomaly]=[a1 0 ];
% 

T1=2*pi()*a1^1.5/mu^0.5;     %Orbital period  

c_w=[cos(argument_of_perigee1) sin(argument_of_perigee1) 0
    -sin(argument_of_perigee1) cos(argument_of_perigee1) 0
    0 0 1];
c_i=[1 0 0
    0 cos(i1) sin(i1)
    0 -sin(i1) cos(i1)];
c_omega=[cos(omega1) sin(omega1) 0
    -sin(omega1) cos(omega1) 0
    0 0 1];
C=c_w*c_i*c_omega;


E=2*atan(tan(true_anomaly1/2)*sqrt((1-e1)/(1+e1)));
r=a1*(1-e1^2)/(1+e1*cos(true_anomaly1));
x=a1*(cos(E)-e1);
y=a1*sqrt(1-e1^2)*sin(E);
xdot=-sqrt(mu*a1)*sin(E)/r;
ydot=sqrt(mu*a1*(1-e1^2))*cos(E)/r;
r_eci1=C'*[x;y;0];
v_eci1=C'*[xdot;ydot;0];
r_eci1=r_eci1';
v_eci1=v_eci1';

% Generate the final orbit position and velocity 
e2=0;
i2=0;
omega2=0;
argument_of_perigee2=90*pi/180;
true_anomaly2=180*pi/180;

% [a,e,i,omega,argument_of_perigee,true_anomaly]=[a1 0 ];
% 

T2=2*pi()*a2^1.5/mu^0.5;     %Orbital period  

c_w=[cos(argument_of_perigee2) sin(argument_of_perigee2) 0
    -sin(argument_of_perigee2) cos(argument_of_perigee2) 0
    0 0 1];
c_i=[1 0 0
    0 cos(i2) sin(i2)
    0 -sin(i2) cos(i2)];
c_omega=[cos(omega2) sin(omega2) 0
    -sin(omega2) cos(omega2) 0
    0 0 1];
C=c_w*c_i*c_omega;


E=2*atan(tan(true_anomaly2/2)*sqrt((1-e2)/(1+e2)));
r=a2*(1-e2^2)/(1+e2*cos(true_anomaly2));
x=a2*(cos(E)-e2);
y=a2*sqrt(1-e2^2)*sin(E);
xdot=-sqrt(mu*a2)*sin(E)/r;
ydot=sqrt(mu*a2*(1-e2^2))*cos(E)/r;
r_eci2=C'*[x;y;0];
v_eci2=C'*[xdot;ydot;0];
r_eci2=r_eci2';
v_eci2=v_eci2';

%% %% Problem III: 1 perform hohmann transfer
[X_initial,X_transfer,X_final]=Hohmann_trans(r_eci1,v_eci1,r_eci2,v_eci2);
X=[X_initial;X_transfer;X_final];

figure(9);
hold on
plot3(X_initial(:,1),X_initial(:,2),X_initial(:,3),'r',LineWidth=2)   %% plot for Problem I-6
hold on
plot3(X_transfer(:,1),X_transfer(:,2),X_transfer(:,3),'g',LineWidth=2)   %% plot for Problem I-6
hold on
plot3(X_final(:,1),X_final(:,2),X_final(:,3),'r',LineWidth=2)   %% plot for Problem I-6
xlabel('ECI x [km]');
ylabel('ECI y [km]');
zlabel('ECI z [km]');
title('Hohmann Transfer using numerical integration');
grid on
hold off

%% Problem III: Perform bi-elliptic transfer
[X_initial,X_transfer,X_final]=bi_elliptic(r_eci1,v_eci1,r_eci2,v_eci2,40*norm(r_eci1));
X=[X_initial;X_transfer;X_final];

figure(10);
hold on
plot3(X_initial(:,1),X_initial(:,2),X_initial(:,3),'r',LineWidth=2)   %% plot for Problem I-6
hold on
plot3(X_transfer(:,1),X_transfer(:,2),X_transfer(:,3),'g',LineWidth=2)   %% plot for Problem I-6
hold on
plot3(X_final(:,1),X_final(:,2),X_final(:,3),'r',LineWidth=2)   %% plot for Problem I-6
xlabel('ECI x [km]');
ylabel('ECI y [km]');
zlabel('ECI z [km]');
title('Hohmann Transfer using numerical integration');
grid on
hold off