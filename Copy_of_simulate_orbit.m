function [positions,velocities] = simulate_orbit(r_eci,v_eci,n)
%SIMULATE_ORBIT simulates a satellite trajectory in cartesian coordinate by solving the kepler's
%equation. 
% INPUT:
%   r_eci => initial position vector in ECI cooridnate system 
%   v_eci => initial velocity vecotr in ECI coordinate system 
%   n     => n is the number of periods over which the simulation is carried out 
%
% OUTPUT: 
%   Positon    => Position[t] is the ECI Position of the satellitecomputed at time t  
%   Velocities => velocities[t] is the ECI velocity of the satellitecomputed at time t  

% Initail values, parameters and variable definitions


mu= 398601.2; 
pi=3.14159265359;

H=cross(r_eci,v_eci);
p=norm(H)^2/mu;
[a,e,i,omega,argument_of_perigee,true_anomaly]=ECI2classical(r_eci,v_eci);    %ECI to classical orbital elements convertion Problem I-1
T=2*pi()*a^1.5/mu^0.5;   %Orbital period  Problem I-2
runtime= n*floor(T);       %for how long to simulate the orbit propagation. Five orbital period length.
t_a=[];                  % t_a[t+1] is the true anomaly at time t
pi_dec=0;


% track true anomaly change for 5 orbital period solving Kepler's equation using Newton Rapson method

% Handles possible issues with quadrant ambiguity
if true_anomaly>=pi
        true_anomaly= true_anomaly-pi;
        pi_dec=1;
end

t=0;
t_a(t+1)=true_anomaly;


while t<runtime 
    E=2*atan(tan(true_anomaly/2)*sqrt((1-e)/(1+e)));
    M= E-e*sin(E);
    delta_t=M*a^1.5/sqrt(mu);
    Tp=t-delta_t;
    t=t+1;
    delta_t=t-Tp;
    M=sqrt(mu/a^3)*delta_t;
    E=solve_Kepler(e,M);
    true_anomaly=2*atan(tan(E/2)*sqrt((1+e)/(1-e)));
    t_a(t+1)=true_anomaly;
end

% applies correction if prior corrections were employed to avoid quadrant
% ambiguity 
if pi_dec
         t_a=t_a+ones(1,length(t_a))*pi();
end




mu= 398601.2; 
pi=3.14159265359;

H=cross(r_eci,v_eci);
p=norm(H)^2/mu;
[a,e,i,omega,argument_of_perigee,true_anomaly]=ECI2classical(r_eci,v_eci);    %ECI to classical orbital elements convertion Problem I-1
T=2*pi()*a^1.5/mu^0.5;   %Orbital period  Problem I-2

t_a=track_true_anomaly(r_eci,v_eci,n);

% simulate the trajectory of satellite in cartesian coordinate 
c_w=[cos(argument_of_perigee) sin(argument_of_perigee) 0
    -sin(argument_of_perigee) cos(argument_of_perigee) 0
    0 0 1];
c_i=[1 0 0
    0 cos(i) sin(i)
    0 -sin(i) cos(i)];
c_omega=[cos(omega) sin(omega) 0
    -sin(omega) cos(omega) 0
    0 0 1];
C=c_w*c_i*c_omega;  

positions=[];
velocities=[];
for i=1:length(t_a)
    E=2*atan(tan(t_a(i)/2)*sqrt((1-e)/(1+e)));
    r=a*(1-e^2)/(1+e*cos(t_a(i)));
%     r_=a*(1-e*cos(E))
    x=a*(cos(E)-e);
%     x_vals(i)=r*cos(t_a(i))
    y=a*sqrt(1-e^2)*sin(E);
%     y_vals(i)=r*sin(t_a(i))
    xdot=-sqrt(mu*a)*sin(E)/r;
%     xdot_vals(i)=-sqrt(mu/p)*sin(t_a(i))
    ydot=sqrt(mu*a*(1-e^2))*cos(E)/r;
%     ydot_vals(i)=sqrt(mu/p)*(e+cos(t_a(i)))
    position=C'*[x;y;0];
    velocity=C'*[xdot;ydot;0];
    positions(i,1:3)=position';
    velocities(i,1:3)=velocity';
end

end

