function g = sqrt(f)
%SQRT   Square root of a CHEBTECH.
%   SQRT(F) returns the square root of a SINGFUN F. Note, it is assumed that F
%   is non-zero on its domain.

% Simply call the compose function:
g = compose(f, @sqrt);

% Throw a warning if the result is not happy and we find roots in the domain:
if ( f.ishappy && ~g.ishappy )
    r = roots(f);
    if ( ~isempty(r) )
        warning(['Attempting to SQRT a CHEBTECH with one or more roots. ', ...
            'Result may be inaccurate.']);
    end
end

end