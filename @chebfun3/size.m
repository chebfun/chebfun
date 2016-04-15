function varargout = size(f , dim)
%SIZE  Size of a CHEBFUN3
%   D = SIZE(F) returns the three-element row vector D = [inf,inf, inf].
%
%   [M, N, P] = SIZE(F) returns M = inf, N = inf and P = inf.
%
%   M = SIZE(F, DIM) returns the dimension specified by the scalar DIM, 
%   which is always inf.

if ( (nargin == 1) && (nargout <= 1) )
    varargout = {[Inf, Inf, Inf]};
    
elseif ( (nargout == 2) && (nargout <=1) )
    varargout = {Inf, Inf, Inf};
    
elseif ( (nargin == 2) && ( (dim == 1) || (dim == 2) || (dim == 3) ) )
    varargout = {Inf};
else
    error('CHEBFUN:CHEBFUN3:size:outputs', 'Too many output arguments.');
end

end