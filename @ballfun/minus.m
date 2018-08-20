function h = minus(f, g)
% MINUS Soustraction between two BALLFUN functions
%   MINUS(f, g) is the soustraction of the BALLFUN function f by g

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

h = f+(-g);
end
