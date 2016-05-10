function e = end(f, k, n)
%END   Rightmost point of a CHEBFUN domain (or last row/col of quasimatrix).
%   END(F, K, N) This function is called for indexing expressions involving a
%   CHEBFUN F when END is part of the K-th index out of N indices.
%
%   If F is a column CHEBFUN and K is 1, this function returns the last point
%   in F's domain. If K is 2, it returns the index of F's last column.
%
%   If F is a row CHEBFUN and K is 1, this function returns the index of F's
%   last row. If K is 2, it returns the last point in F's domain.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( n > 2 )
    error('CHEBFUN:CHEBFUN:end:ngt2', 'Index exceeds CHEBFUN dimensions.');
end

if ( ((k == 2) && ~f(1).isTransposed) || ...
     ((k == 1) && f(1).isTransposed && (n > 1)) )
    % 'end' row/column of the array-valued CHEBFUN.
    if ( isempty(f) )
        e = 0;
    else
        e = numColumns(f);
    end
    
else
    % 'end' of the domain.
    e = f(1).domain(end);
    
end

end
