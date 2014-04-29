function t = isFactored(disc)
% ISFACTORED
% 
% T = ISFACTORED(DISC) returns true if the discretization has stored LU factors
% whose dimensions are compatible with the current discretization size.

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Assume false
t = false;

% Check whether any factorisation data is stored
if ( ~isempty(disc.mldivideData) )
    L = disc.mldivideData{1};
    % Check whether dimensions stored are compatible with current discretization
    % size.
    if ( sum(disc.dimension) == size(L, 1) )
        t = true;
    end
end

end
