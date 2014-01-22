function out = sum(f)
%SUM   Definite integral of a DELTAFUN.
%   SUM(F) is the integral of F on the domain of F.
%
% See also CUMSUM, DIFF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.


%%
% Trivial case:
if ( isempty(f) )
    % Following Matlab:
    out = 0;
    return
end

out = sum(f.funPart);

% Add integral of delta functions:
if ( ~isempty(f.deltaMag) )
    % What happens for the higher order derivatives etc?
    % Answer: Neglect higher order deltas, since the integral of any derivative
    % of a delta function can be considered as it's action on the function 1 and
    % hence it's integral is zero.
    deltaMag = f.deltaMag;
    out = out + sum(deltaMag(1, :));
end    
end