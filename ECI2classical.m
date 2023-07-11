function [a,e,i,omega,argument_of_perigee,true_anomaly] = ECI2classical(position,velocity)
%this function converts position and velocity in ECI coordinate to for a
%given eliptical orbit. Position and velocity must be given as row vector. 
%classical orbital elements [a,e,i,omega,argument_of_perigee,true_anomaly]

mu= 398601.2; 
pi=3.14159265359;
r=sqrt(position(1)^2 +position(2)^2+ position(3)^2);   % scalar norm of positon 
v=sqrt(velocity(1)^2 + velocity(2)^2+ velocity(3)^2);  % scalar norm of velocity 

H=cross(position,velocity);    % vector H
h= norm(H);           %scalar norm of H

a= 1/(-v^2/mu + 2/r);
p= h^2/mu;

e_vec= 1/mu*(cross(velocity,H)-mu*position/r);  % eccentricity vector 
e= norm(e_vec);  %eccentricity 

u_H= H/h;  %unit vector in H direction 
u_e= e_vec/e; %unit vector in e direction 
u_m=cross(u_H,u_e);

i = acos(u_H(3));
% argument_of_perigee= atan(u_H(1)/u_m(3));


n= cross([0 0 1],H);
n= n/norm(n);

% if dot(n,e_vec)<0 
%     argument_of_perigee= 2*pi()-acos(dot(u_e,position/norm(position)));
% else 
%     argument_of_perigee=acos(dot(u_e,position/norm(position)));
% end


% if dot(position,velocity)>0
%     true_anomaly=acos(dot(u_e,position/norm(position)));
%     if ~(true_anomaly<pi()) & true_anomaly<0 
%             true_anomaly=true_anomaly+2*pi;
% 
%     elseif ~(true_anomaly<pi()) & true_anomaly>0 
%             true_anomaly=-true_anomaly+2*pi;
%     end
% else 
%    true_anomaly=acos(dot(u_e,position/norm(position)));
%     if ~(true_anomaly>pi()) & true_anomaly<0 
%             true_anomaly=true_anomaly+2*pi;
%     elseif ~(true_anomaly>pi()) & true_anomaly>0 
%             true_anomaly=-true_anomaly+2*pi;
%     end

if dot(position,velocity)>0
    true_anomaly=acos(dot(u_e,position/norm(position)));
else 
   true_anomaly=2*pi-acos(dot(u_e,position/norm(position)));
end


% argument_of_perigee= atan(u_H(1)/u_m(3))*180/pi

% if  u_H(1)<0 & u_m(3)<0
%     argument_of_perigee= atan(u_H(1)/u_m(3))+pi;
% elseif u_H(1)<0 & u_m(3)>0;
%     argument_of_perigee= atan(u_H(1)/u_m(3))+pi;
% elseif u_m(3)<0 & u_H(1)>0;
%     argument_of_perigee= atan(u_H(1)/u_m(3))+2*pi;
% else 
%     argument_of_perigee= atan(u_H(1)/u_m(3));
% end



% 

argument_of_perigee= acos(dot(n,u_e));
if e_vec(3)>0
    argument_of_perigee= acos(dot(n,u_e));
else 
    if acos(dot(n,u_e))>0   
         argument_of_perigee= 2*pi()-acos(dot(n,u_e));
    end
end



if  n(1)<0 & n(2)<0
    omega= atan(n(2)/n(1))+pi;
elseif n(1)<0 & n(2)>0;
    omega= atan(n(2)/n(1))+pi;
elseif n(2)<0 & n(1)>0;
    omega= atan(n(2)/n(1))+2*pi;
else 
    omega= atan(n(2)/n(1));
end


% classical= [a,e,i*180/pi,omega*180/pi,argument_of_perigee*180/pi,true_anomaly*180/pi]

end





