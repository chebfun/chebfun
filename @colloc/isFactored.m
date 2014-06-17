function t = isFactored(disc)
%ISFACTORED   Check if a useful factorization of a COLLOC is stored.
%   ISFACTORED(DISC) returns true if the discretization has stored LU
%   factors whose dimensions are compatible with the current discretization
%   size.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Assume false:
t = false;

% Check whether any factorisation data is stored:
if ( ~isempty(disc.mldivideData) )
    % Total discretization size is the sum over the pieces, times the
    % number of function variables, plus the number of scalar variables. 
    isFun = isFunVariable(disc.source);
    len = sum(isFun)*sum(disc.dimension) + sum(~isFun);
    
    % Compare length to that of stored data.
    L = disc.mldivideData{1};
    if ( len == size(L, 1) )
        t = true;
    end
end

end
