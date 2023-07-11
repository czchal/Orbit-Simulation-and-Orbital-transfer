%Initial orbit

%circular orbit of altitude 200km

r_eci1=[1466.16169146486   -6305.61776724743    -1166.04180684668];
v_eci1=[4.49930127355746   2.15141507694868     -5.97692840415413];

%final circular target orbit of geo

r_eci2=[9400 -40422 -7477];
v_eci2=[1.777074590403157   0.8497375112720199  -2.360688237994593];


param.mu = 398601.2;          %[km^3/s^2]
param.Re = 6378;  

a1=param.Re+200;

v1=sqrt(param.mu/a1);


a2=42167;

v2=sqrt(param.mu/a2);

%transfer orbit semi major axis 
at=(a1+a2)/2;

%first burn
deltav1=sqrt(2*param.mu*(1/a1-1/(a1+a2)))-sqrt(param.mu/a1)
%second burn
deltav2=-sqrt(2*param.mu*(1/a2-1/(a1+a2)))+sqrt(param.mu/a2)


%simulate the initial orbit for one period 

X_initial=num_orbit_nodrag(r_eci1,v_eci1,1,param);



% after the first burn the satellite travels half period of the ellipse
% until apogee


v=norm(v_eci1);   % 

% Simulates Hohmann transfer between two coplanar orbits 

v=v+deltav1;

v_vec=v*v_eci1/norm(v_eci1);

v=norm(v_vec);
r=norm(r_eci1);

a= 1/(-v^2/param.mu + 2/r);

T=2*pi()*a^1.5/param.mu^0.5;

X_transfer=num_orbit_nodrag(r_eci1,v_vec,0.5,param);



%simulate the final orbit for one period




v_final=X_transfer(end,4:6)

v_c2=norm(v_final)+deltav2;
v_c2_vec=v_c2*v_final/norm(v_final)

r_e=X_transfer(end,1:3)   
X_final=num_orbit_nodrag(r_e,v_c2_vec,1,param);

X=[X_initial;X_transfer;X_final];

figure(1);
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
