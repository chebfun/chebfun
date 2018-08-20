function w = power(v,n)
% POWER Return a BALLFUNV to the power n
%   POWER(v) is the BALLFUNV (Vx^n, Vy^n, Vz^n)

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

V = v.comp;
w = ballfunv(power(V{1},n),power(V{2},n),power(V{3},n));
end
