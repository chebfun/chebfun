function varargout = x
%X   A chebfun of the identity on [-1,1].
%   X = CHEB.X returns a chebfun object for the function @(x)x on [-1,1].
%
%   CHEB.X is shorthand for the expression X = CHEBFUN(@(X) X).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

x = chebfun(@(x) x);

if ( nargout > 0 ) 
    
    % For the syntax x = cheb.x:
    varargout{1} = x; 
    
    if ( nargout > 1 ) 
        error('CHEB:X:TooManyOutputs',... 
            'Too many output arguments. CHEB.X only returns "x".')
    end
    
else
    
   % Put 'x' into the workspace:
   assignin('base', 'x', x)
   
end
end
