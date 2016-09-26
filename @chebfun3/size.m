function varargout = size(f, dim)
%SIZE   Size of a CHEBFUN3 object.
%   D = SIZE(F) returns the three-element row vector D = [Inf, Inf, Inf].
%
%   [M, N, P] = SIZE(F) returns M = Inf, N = Inf and P = Inf.
%
%   M = SIZE(F, DIM) returns the size of F along the dimension specified by
%   the scalar DIM, and is always Inf.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( (nargin == 1) && (nargout <= 1) )
    varargout = {[Inf, Inf, Inf]};
    
elseif ( (nargin <=1) && (nargout == 3) )
    varargout = {Inf, Inf, Inf};
    
elseif ( (nargin == 2) && ( (dim == 1) || (dim == 2) || (dim == 3) ) )
    varargout = {Inf};
else
    error('CHEBFUN:CHEBFUN3:size:outputs', 'Too many output arguments.');
end

end