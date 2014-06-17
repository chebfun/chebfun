function out = sum(f)
%SUM   Definite integral of a DELTAFUN.
%   SUM(F) is the integral of F over the domain of F.
%
% See also CUMSUM, DIFF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = sum(f.funPart);

% Add integral of delta functions:
if ( ~isempty(f.deltaMag) )
    % What happens for the higher order derivatives etc? Answer: Neglect higher
    % order deltas since the integral of any derivative of a delta function can
    % be considered as it's action on the function 1. Hence, it's integral is 0.
    out = out + sum(f.deltaMag(1, :));
end    

end
