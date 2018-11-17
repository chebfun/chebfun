function w = real(v)
% REAL  Real part of a BALLFUNV
%   REAL(V) is the real part of the BALLFUNV V.
%
% See also IMAG.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

V = v.comp;
w = ballfunv( real(V{1}), real(V{2}), real(V{3}) );
end
