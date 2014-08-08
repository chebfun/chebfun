function [A, B] = domain(f, flag)
%DOMAIN   Domain of definition of a CHEBFUN.
%   I = DOMAIN(F) returns a row vector representing the domain of definition
%   (including breakpoints) of the CHEBFUN F, and is equivalent to F.domain.
% 
%   [A, B] = DOMAIN(F) returns the endpoints of the domain as scalars and I =
%   DOMAIN(F, 'ENDS') returns a vector of the end points.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargout == 2 )
    % Return the end points as two scalars:
    dom = f(1).domain([1 end]);
    A = dom(1);
    B = dom(2);

elseif ( nargin == 2 )
    % Return the end points as a vector
    if ( strcmpi(flag, 'ends') )
        A = f(1).domain([1 end]);
    else
        error('CHEBFUN:CHEBFUN:domain:unknown', 'Unexpected input.');
    end
    
elseif ( numel(f) == 1 )
    % CHEBFUN case:
    A = f.domain;
    
else
    % Merge the domains of columns in a quasimatrix:
    dom = cell(1, numel(f));
    for k = 1:numel(f)
        dom{k} = f(k).domain;
    end
    A = domain.merge(dom{:});
end

end
