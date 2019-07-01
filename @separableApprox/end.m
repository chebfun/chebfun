function e = end(f, k, n)
%END   Rightmost point of a separableApprox domain.
%   END(F, K, N) This function is called for indexing expressions involving a
%   separableApprox F when END is part of the K-th index out of N indices.
%
%   If K is 1, this function returns the last point in F's domain in the 'x'
%   direction. If K is 2, it returns the index of F's domain in the 'y'
%   direction

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( n > 2 )
    error('CHEBFUN:separableApprox:end:ngt2', ...
        'Index exceeds separableApprox dimensions.');
end

dom = f.domain;

if ( k == 1 )
    
    % 'end' of the x domain.
    e = dom(2);
    
elseif ( k == 2 )
    
    % 'end' of the y domain.
    e = dom(4);
    
end

end
