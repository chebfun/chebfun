function w = imag(v)
%IMAG   Complex imaginary part of a BALLFUNV
%   IMAG(V) is the imaginary part of the BALLFUNV V.
%
% See also REAL.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty( v )
    w = ballfunv();
    return
end

V = v.comp;
w = [ imag(V{1}) ; imag(V{2}) ; imag(V{3})];
end
