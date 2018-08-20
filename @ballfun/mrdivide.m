function g = mrdivide(f,c)
% MRDIVIDE / Right scalar divide for BALLFUN object
%   MRDIVIDE(f,c) divides the BALLFUN function f by the scalar c

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = f*(1/c);
end
