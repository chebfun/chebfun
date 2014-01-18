function f = fred( K, v )
%FRED  Fredholm integral operator with a chebfun2 kernel.
%
%  F = FRED(K, V) computes the Fredholm integral with kernel K:
%
%       (F*v)(x) = int( K(x,y)*v(y), y=c..d ),  x=a..b
%
%  where [c d] = domain(V) and [a b c d] = domain(K). The kernel function
%  K(x,y) should be smooth for best results. K is a chebfun2 and V is a chebfun.
%  The result is a row chebfun object.
%
% See also VOLT.

% Copyright 2013 by The University of Oxford and The Chebfun2 Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check:
if ( isempty( K ) )
    return
end

% Get domain of f:
dom = K.domain;

% Get the low rank representation for f: 
cols = K.cols; 
rows = K.rows; 
piv = K.pivotValues; 
d = 1./piv; 
d(d==inf) = 0;  % set infinite values to zero. 

% Make sure v is a chebfun:
if ( isa( v, 'function_handle' ) )
    v = chebfun( v, dom(3:4) );
else
    % Check that the domain of chebfun2 and chebfun are correct.
    if ( ~domainCheck( v, cols ) )
        error('CHEBFUN2:FRED:CHEBDOMAIN','Domain of chebfun and chebfun2 kernel do not match');
    end
end

cols = ( cols.' * v ).';
f = cols * diag( d ) * rows.';

end