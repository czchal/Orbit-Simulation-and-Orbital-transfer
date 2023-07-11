function dx = twobody_drag(t, x, param)
% ad= 0.5*rho*cd*S*v^2/m
r = x(1:3);
v = x(4:6);

u_v= v/norm(v);   %velocity direction
ad= 0.5*param.rho*param.cd*param.S*norm(v)^2/param.m;  % drag due to air 

dr = v;
dv = - param.mu / norm(r)^3 * r - ad*u_v;

dx = [dr, dv];

end


% exp(-(norm(X(1,1:3))-param.Re)/7)