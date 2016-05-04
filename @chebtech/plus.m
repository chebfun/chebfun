function f = plus(f, g)
%+   Addition of two CHEBTECH objects.
%   F + G adds F and G, where F and G may be CHEBTECH objects or scalars.
%
%   If F is an array-valued CHEBTECH, then F + C is supported if C is a row
%   vector of doubles with the same number of columns as F.
%
% See also MINUS, UPLUS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(g) ) % CHEBTECH + [] = []
    
    f = [];
    
elseif ( isa(g, 'double') ) % CHEBTECH + double
    % Update values (use bsxfun() to handle the case in which g is a vector
    % and f is an array-valued CHEBTECH):
    % Update coeffs:
    if ( (size(g, 2) > 1) && (size(f.coeffs, 2) == 1) )
        % Perform singleton expansion of f:
        f.coeffs = repmat(f.coeffs, 1, size(g, 2));
    end
    f.coeffs(1,:) = f.coeffs(1,:) + g;
elseif ( isa(f, 'double') ) % double + CHEBTECH
    
    % Switch argument order and call CHEBTECH/PLUS again:
    f = plus(g, f);
    
elseif ( isa(f, 'chebtech') && isa(g, 'chebtech') )  % CHEBTECH + CHEBTECH
    
    % Make both CHEBTECH objects have the same length:
    nf = size(f.coeffs, 1);
    ng = size(g.coeffs, 1);
    if ( nf > ng )
        % Increase the length of g (via PROLONG):
        g = prolong(g, nf);
    elseif ( nf < ng )
        % Increase the length of f (via PROLONG):
        f = prolong(f, ng);
    end

    % Tolerance for determining zero-output:
    tol = eps*max(vscale(f), vscale(g));

    % Update values and coefficients:
    f.coeffs = f.coeffs + g.coeffs;
    
    % Look for a zero output:
    absCoeffs = abs(f.coeffs);
    isz = bsxfun(@lt, absCoeffs, .2*tol); % Are coeffs below .2*el*vs?
    
    if ( all(isz(:)) )
        % Create a zero CHEBTECH:
        ishappy = f.ishappy && g.ishappy;
        z = zeros(1, size(f.coeffs, 2));

        data.vscale = z;
        f = f.make(z, data);
        f.ishappy = ishappy;
    else
        f.ishappy = f.ishappy && g.ishappy;
    end

else    % Don't know how to do the addition of the objects
    
    error('CHEBFUN:CHEBTECH:plus:typeMismatch', ...
        ['Incompatible operation between objects.\n', ...
         'Make sure functions are of the same type.']);
    
end

end
