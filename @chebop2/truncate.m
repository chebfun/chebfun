function a = truncate(a, tol)
%TRUNCATE   Truncation of CHEBOP2 objects.
%   A = TRUNCATE(A, TOL) truncates off all elements in A which are below absolute
% value TOL.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Fast binary search.
factor = 2;

while ( 1 )
    cut = floor(length(a)*(1-1/factor));
    if ( cut < 1 || length(a) <= 1 )
        break
    end
    if ( max(abs(a(cut:end))) < tol )
        a = a(1:cut);
    else
        factor = factor*2;
    end
    if ( factor > 16 )
        break
    end
end

if ( length(a) <= 1 )
    return
end

% Slower truncation to finish the job.
while ( abs(a(end)) < tol )
    a(end) = [];
end

end