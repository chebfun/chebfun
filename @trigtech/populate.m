function f = populate(f, op, data, pref)
%POPULATE   Populate a TRIGTECH class with values.
%   F = F.POPULATE(OP) returns a TRIGTECH representation populated with values
%   F.VALUES of the function OP evaluated on an equally spaced grid. The field
%   F.ISHAPPY indicates whether the representation is deemed 'happy'.
%   Essentially this means that such an interpolant is a sufficiently accurate
%   approximation to OP. If F.ISHAPPY is FALSE, then POPULATE was not able to
%   obtain a happy result.
%
%   OP should be vectorized (i.e., accept a vector input), and output a vector
%   of the same length. Furthermore, OP may be an array-valued function, in
%   which case it should accept a vector of length N and return a matrix of size
%   NxM.
%
%   F.POPULATE(OP, VSCALE, HSCALE) enforces that the happiness check is relative
%   to the initial vertical scale VSCALE and horizontal scale HSCALE. These
%   values default to 0 and 1 respectively. During refinement, VSCALE updates
%   itself to be the largest magnitude values to which (each of the columns in)
%   OP evaluated to.
%
%   F.POPULATE(OP, VSCALE, HSCALE, PREF) enforces any additional preferences
%   specified in the preference structure PREF (see TRIGTECH.TECHPREF).
%
%   F.POPULATE(VALUES, ...) (or F.POPULATE({VALUES, COEFFS}, ...)) populates F
%   non-adaptively with the VALUES (and COEFFS) passed. These values are still
%   tested for happiness in the same way as described above, but the length of
%   the representation is not altered.
%
% See also TRIGTECH, TECHPREF, HAPPINESSCHECK.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%% Non-adaptive construction. %%%%%%%%%%%%%%%%%%%%%%%%%%
% Values (and possibly coefficients) have been given.
if ( isnumeric(op) || iscell(op) )
    if ( isnumeric(op) )
        % OP is just the values.
        f.values = op;
        f.coeffs = f.vals2coeffs(op);
    else                 
        % OP is a cell {values, coeffs}
        f.coeffs = op{2};
        f.values = f.coeffs2vals(f.coeffs);
    end

    % We're always happy if given discrete data:
    f.ishappy = true;

    % Threshold to determine if f is real or not.
    vscl = max(data.vscale, vscale(f));
    f.isReal = max(abs(imag(f.values))) <= 3*(eps*vscl);
    f.values(:,f.isReal) = real(f.values(:,f.isReal));

    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Adaptive construction. %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise empty values to pass to refine:
f.values = [];
% Check a random value of op in (-1,1) to see if the result is complex and
% if the value is a NaN or Inf.
pseudoRand = 0.376989633393435;
rndVal = feval(op, (2*pseudoRand - 1));

if ( any(isnan(rndVal)) || any(isinf(rndVal)) )
    error('CHEBFUN:TRIGTECH:populate:isNan','Cannot handle functions that evaluate to Inf or NaN.');
end

f.isReal = false(size(rndVal));
for k = 1:numel(rndVal)
    f.isReal(k) = isreal(rndVal(k));
end

% Loop until ISHAPPY or GIVEUP:
while ( 1 )

    % Call the appropriate refinement routine: (in PREF.REFINEMENTFUNCTION)
    [f.values, giveUp] = f.refine(op, f.values, pref);

    % We're giving up! :(
    if ( giveUp ) 
        break
    end    
    
    % Update vertical scale: (Only include sampled finite values)
    valuesTemp = f.values;
    valuesTemp(~isfinite(f.values)) = 0;
    data.vscale = max(data.vscale, max(abs(valuesTemp)));
    
    % Compute the trigonometric coefficients:
    coeffs = f.vals2coeffs(f.values);
    
    % Check for happiness:
    f.coeffs = coeffs;
    [ishappy, cutoff] = happinessCheck(f, op, f.values, data, pref);
    
        
    % We're happy! :)
    if ( ishappy ) 
        f = prolong(f,cutoff); % chop coefficients to length cutoff
        break
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%% Assign to TRIGTECH object. %%%%%%%%%%%%%%%%%%%%%%%%%%
f.ishappy = ishappy;

% Force the values to be real if the imaginary part is zero
f.values(:,f.isReal) = real(f.values(:,f.isReal));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ouput. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( ishappy )
    % We're done, and can return.
    return
end

end
