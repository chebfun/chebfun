function f = cheb2four(f)
%CHEB2FOUR   Convert a Chebyshev-based CHEBFUN to a Fourier-based CHEBFUN.
%   G = CHEB2FOUR(F) converts the Chebyshev-based CHEBFUN F to a Fourier-based
%   CHEBFUN G.
%
% See also FOUR2CHEB.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% TODO: This violates encapsulation in a big way..

if ( isempty(f) )
    f = chebfun([], 'periodic');
    
elseif ( isa(f.funs{1}.onefun, 'fourtech') )
    % Nothing to do.

elseif ( isa(f.funs{1}.onefun, 'chebtech') )
    f = chebfun(f,f.domain, 'periodic');
    
else
    error('CHEBFUN:cheb2four:TypeUnkown', 'Tech type unknown');
    
end

end