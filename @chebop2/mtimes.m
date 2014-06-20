function N = mtimes(N, u)
%MTIMES   Forward application of a CHEBOP2 object.
% 
% N = mtimes(N, u) is the same as N evaluated at the CHEBFUN2 u. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(u, 'chebfun2') )
    op = N.op;
    N = op(u);     
elseif ( isa(N, 'chebop2') && isa(u, 'double') )
    N.coeffs = u*N.coeffs; 
    op = N.op; 
    N.op = @(v) u*op(v);    
elseif ( isa(u, 'chebop2') && isa(N, 'double') )
    N = mtimes(u, N);
else
    error('CHEBFUN:CHEBOP2:mtimes:badInput', ...
        'Can only times a CHEBOP2 by a double or forward apply to a CHEBFUN2.');
end

end
