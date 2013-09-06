function e = end(f, k, n)
%END   Rightmost point of a CHEBFUN domain (or last row/col of quasimatrix).

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( n > 2 )
    error('CHEBFUN:end:ngt2', 'Index exceeds CHEBFUN dimensions.');
end


if ( k == 2 && ~f.isTransposed ) || ( k == 1 && f.isTransposed )
    % 'end' row/column of the quasimatrix.
    if ( isempty(f) )
        e = 0;
    else
        e = size(f.funs, 2);
    end
    
else
    % 'end' of the domain.
    e = f.domain(end);
    
end

end