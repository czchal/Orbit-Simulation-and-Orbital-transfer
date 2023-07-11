function [X_initial,X_transfer,X_final] = Hohmann_trans(r_eci1,v_eci1,r_eci2,v_eci2)
%Hohmann Trans computes the trajectory of hohmann transfer for two circular
%orbits using RK-4 numerical integration
%  
% Inputs:
%   r_eci1 = the eci position coordinate of the spacraft at the inital transfer point
%   v_eci1 = the eci velocity coordinate of the spacraft at the inital transfer point
%   r_eci2 = the eci position coordinate of the spacraft at the the transfer destination 
%   v_eci2 = the eci velocity coordinate of the spacraft at the the transfer destination
%
% Output:
%   X_initial = Trajectory of the inital orbit for one Period
%   X_transfer=  Trajectory of the transfer orbit 
%   X_final = Trajectory of the final orbit for one Period

param.mu = 398601.2;          %[km^3/s^2]
param.Re = 6378;  

a1=norm(r_eci1)

v1=sqrt(param.mu/a1);


a2=norm(r_eci2);

v2=sqrt(param.mu/a2);

%transfer orbit semi major axis 
at=(a1+a2)/2;

%first burn
deltav1=sqrt(2*param.mu*(1/a1-1/(a1+a2)))-sqrt(param.mu/a1)
%second burn
deltav2=-sqrt(2*param.mu*(1/a2-1/(a1+a2)))+sqrt(param.mu/a2)


%simulate the initial orbit for one period 

X_initial=num_orbit_nodrag(r_eci1,v_eci1,1,param);     % trajectory of the inital orbit for one Period



% after the first burn the satellite travels half period of the ellipse
% until apogee
% Simulates Hohmann transfer between two coplanar orbits 



v=norm(v_eci1);       % magnitude of velocity at initial orbit 

v=v+deltav1;          % magnitude of velocity of the transfer orbit at the perigee point

v_vec=v*v_eci1/norm(v_eci1);

v=norm(v_vec);
r=norm(r_eci1);

a= 1/(-v^2/param.mu + 2/r);     %semi major axis of the transfer orbit

T=2*pi()*at^1.5/param.mu^0.5;    %period of the transfer orbit


X_transfer=num_orbit_nodrag(r_eci1,v_vec,0.5,param);    % Trajectory of the transfer orbit



%simulate the final orbit for one period



v_final=X_transfer(end,4:6)   % The velocity of the transfer orbit at apogee 

v_c2=norm(v_final)+deltav2;    % magnitude of Velocity on the final orbit

v_c2_vec=v_c2*v_final/norm(v_final)
 
r_e=X_transfer(end,1:3)          

X_final=num_orbit_nodrag(r_e,v_c2_vec,1,param);

% X=[X_initial;X_transfer;X_final];


end

