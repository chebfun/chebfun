function w = real(v)
%REAL  Real part of a BALLFUNV.
%   REAL(F) returns the BALLFUNV representing the real part of F.
%
% See also IMAG.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty( v )
    w = ballfunv();
    return
end

% Take real part of each component:
w = ballfunv( real(v.comp{1}), real(v.comp{2}), real(v.comp{3}) );
end
