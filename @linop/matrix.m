function A = matrix(L,varargin)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
dsc = L.discretizer( L, varargin{:} );
A = matrix(dsc);
end
