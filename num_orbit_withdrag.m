function [X,semi_major] = num_orbit_withdrag(r_eci,v_eci,n,param)
%NUM_ORBIT_NODRAG simulates orbit by numerically integrating the Govering
%equation considering the aerodynamic drag on the spacecraft 
pi=3.14159265359;

%%for inital condition of r and v calculate the required parameters like a, T and dt

r=norm(r_eci);
v=norm(v_eci);
a= 1/(-v^2/param.mu + 2/r);
T=2*pi()*a^1.5/param.mu^0.5;

dt=1;

t = 0:dt:n*floor(T);               
n = length(t);

X = zeros(n, 6);
X(1, :) = [r_eci v_eci];


h = dt;
semi_major=a;
for i = 1:1:n-1
    k1 = h * twobody_drag(t(i), X(i, :), param);
    k2 = h * twobody_drag(t(i) + 0.5 * h, X(i, :) + 0.5 * k1, param);
    k3 = h * twobody_drag(t(i) + 0.5 * h, X(i, :) + 0.5 * k2, param);
    k4 = h * twobody_drag(t(i) + h, X(i, :) + k3, param);
    X(i+1, :) = X(i, :) + 1 / 6 * (k1 + 2*k2 + 2*k3 + k4);
    r=norm(X(i+1,1:3));
    v=norm(X(i+1,4:6));
    a= 1/(-v^2/param.mu + 2/r);
    semi_major(i+1)=a;
end

