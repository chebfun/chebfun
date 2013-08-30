function [A, B] = domain(f)
%DOMAIN   Domain of definition of a CHEBFUN.
%   I = DOMAIN(F) returns the a row vector representing the domain of definition
%   (including breakpoints) of the CHEBFUN F, and is equivalent to F.domain.
% 
%   [A, B] = DOMAIN(F) returns the endpoints of the domain as scalars.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargout == 1)     % One output.
    A = f.domain;
elseif ( nargin == 2 ) % Two outputs.
    A = f.domain(1);
    B = f.domain(end);
end

end
