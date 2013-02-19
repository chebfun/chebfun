function [siz1, siz2] = size(f, varargin)
%SIZE	Size of a FUNCHEB1.
%   SIZE(F) is the number of values at Chebyshev points used to define F.
%
%   See also LENGTH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% The size of a FUNCHEB2 is the size of its vector of values.
siz1 = size(f.values, varargin{:});

if ( nargout == 2 )
    siz2 = siz1(2);
    siz1 = siz1(1);
end

end
