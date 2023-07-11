function dx = twobody(t, x, param)

r = x(1:3);
v = x(4:6);

dr = v;
dv = - param.mu / norm(r)^3 * r;

dx = [dr, dv];

end

