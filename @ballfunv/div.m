function f = div(v)
% DIV Divergence of a BALLFUNV in cartesian system
%   DIV(f) is the divergence of the BALLFUNV v expressed in the
%   cartesian system

[Vx,Vy,Vz] = v.comp{:};
f = diff(Vx,1,"cart") + diff(Vy,2,"cart") + diff(Vz,3,"cart");
end