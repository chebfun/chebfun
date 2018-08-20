function h = cross(f,g)
% CROSS Cross product of two BALLFUNV
%   CROSS(f, g) is the cross product of the BALLFUNV f and g

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.comp;
G = g.comp;
% Get the values of f and g
F1 = ballfun.coeffs2vals(F{1}.coeffs);
F2 = ballfun.coeffs2vals(F{2}.coeffs);
F3 = ballfun.coeffs2vals(F{3}.coeffs);
G1 = ballfun.coeffs2vals(G{1}.coeffs);
G2 = ballfun.coeffs2vals(G{2}.coeffs);
G3 = ballfun.coeffs2vals(G{3}.coeffs);

% Multiply the point values
h1 = ballfun.vals2coeffs(F2.*G3 - F3.*G2);
h2 = ballfun.vals2coeffs(F3.*G1 - F1.*G3);
h3 = ballfun.vals2coeffs(F1.*G2 - F2.*G1);

% Create the ballfunv
h = ballfunv(ballfun(h1),ballfun(h2),ballfun(h3));
end