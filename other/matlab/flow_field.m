function v = flow_field(x, xc, u, w, a)
% This function calculates the flow at the point x in the fluid
% due to a sphere of radius a centred at xc moving with velocity u and
% angular velocity w. This result only holds for an isolated sphere, but it
% will be used as an approximation of the controbution a given sphere makes
% to the overall flow field at a point.

rvec = x - xc;
r = norm(rvec);

% translation
v = (1 + a*a/(3*r*r))*u/r;
v = v + (1 - a*a/(r*r))*(rvec' * u)*rvec/(r*r*r);
v = 3*a*v/4;

% rotation
v = v + a*a*a*cross(w, rvec)/(r*r*r);

end

