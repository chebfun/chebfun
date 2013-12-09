function [A, B] = domain(f)
%DOMAIN   Domain of definition of a CHEBFUN.
%   I = DOMAIN(F) returns a row vector representing the domain of definition
%   (including breakpoints) of the CHEBFUN F, and is equivalent to F.domain.
% 
%   [A, B] = DOMAIN(F) returns the endpoints of the domain as scalars.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Merge the domains of columns in a quasimatrix:
dom = cell(1, numel(f));
for k = 1:numel(f)
    dom{k} = f(k).domain;
end
dom = chebfun.mergeDomains(dom{:});

% Format the output:
if ( nargout <= 1)      % One output.
    A = dom;
elseif ( nargout == 2 ) % Two outputs.
    A = dom(1);
    B = dom(end);
end


end
