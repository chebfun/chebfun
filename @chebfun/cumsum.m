function f = cumsum(f, m)
%CUMSUM   Indefinite integral of a CHEBFUN.
%   G = CUMSUM(F) is the indefinite integral of the CHEBFUN F. G will typically
%   be normalised so that G(F.domain(1)) = 0.
%
%   CUMSUM(F, N) returns the Nth integral of F. If N is not an integer then
%   CUMSUM(F, N) returns the fractional integral of order N as defined by the
%   Riemann-Liouville integral.
%
% See also SUM, INTEGRAL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% [TODO]: Update the above help text once we have singfun.
%   G will typically be normalised so that G(F.domain(1)) = 0.  The exception to
%   this is when computing indefinite integrals of functions which are not
%   integrable at the left boundary. In this case, the arbitrary constant in the
%   indefinite integral is chosen to make the representation of G as simple as
%   possible. Dirac deltas already existing in F will decrease their degree.

% Trivial case:
if ( isempty(f) )
    return
end

% Parse inputs:
if ( nargin == 1 )
    m = 1;
end
    
if ( round(m) ~= m )
    % Fractional integral:
    % [TODO]: Implement this!
    error('CHEBFUN:cumsum:notImplemented', ...
        'Fractional antiderivatives not yet implemented.');
    f = fracCalc(f, m);
    return
end

% Get some basic information from f:
dom = f.domain;
funs = f.funs;
numFuns = numel(funs);
numCols = size(f.funs{1}, 2);

% Loop m times:
for l = 1:m
    
    % Get the level 2 (delta function) impulse data:
    if ( size(f.impulses, 3) > 1 )
        deltas = f.impulses(:,:,2);
    else
        deltas = zeros(length(dom), numCols);
    end
    
    fa = deltas(1,:);
    for j = 1:numFuns
        
%         % [TODO]: Replace this when SINGFUN is added.
%         cumsumFunJ = cumsum(funs{j});
%         if ( nFuns > 1 )
%         % This is because unbounded functions may not be zero at left.
%             lval = get(cumsumFunJ, 'lval');
%             if ( ~isinf(lval) && ~isnan(lval) )
%                 cumsumFunJ = cumsumFunJ - lval;
%             end
%         end
%         funs{j} = cumsumFunJ + fa;
        
        funs{j} = cumsum(funs{j}) + fa;
        fa = get(funs{j}, 'rval') + deltas(j+1,:);
    end
    
    % Get the new impulse data:
    newImps = chebfun.getValuesAtBreakpoints(funs, dom);
    f.impulses = cat(3, newImps, f.impulses(:,:,3:end));
    
end

% Append the updated FUNs:
f.funs = funs;

end

