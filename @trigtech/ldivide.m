function f = ldivide(f, c, varargin)
%.\   Left array divide for a TRIGTECH.
%   C .\ F is equivalent to F ./ C.
%
% See also MLDIVIDE, RDIVIDE, TIMES.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = rdivide(c, f, varargin{:});

end
