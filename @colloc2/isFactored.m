function t = isFactored(disc)

% Returns true if the discretization has stored LU factors whose dimensions are
% compatible with the current discretization size. 

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

t = false;
if ( ~isempty(disc.mldivideData) )
    L = disc.mldivideData{1};
    if ( sum(disc.dimension) == size(L,1) )
        t = true;
    end
end

end
