function v = conj(v)
%CONJ Complex conjugate of a BALLFUNV.
%   CONJ(F) returns the complex conjugate of F. For a complex F, CONJ(F) =
%   REAL(F) - i*IMAG(F).
%
%   See also IMAG, REAL.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty( v )
    return
end

% Take conjugate part of each component:
v = ballfunv( conj(v.comp{1}), conj(v.comp{2}), conj(v.comp{3}) );
end
