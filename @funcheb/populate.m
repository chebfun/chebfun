function f = populate(f, op, vscale, hscale, pref)
%POPULATE    Populate a FUNCHEB class with values.
%   F = F.POPULATE(OP) returns a FUNCHEB representation populated with values
%   F.values of the function OP evaluated on a Chebyshev grid. The fields
%   F.ishappy and F.epslevel denote whether the representation is deemed 'happy'
%   (see HAPPINESSCHECK.m). Essentially this means that such an interpolant is a
%   sufficiently accurate (i.e., to a relative accuracy of F.epsvlevel)
%   approximation to OP. If F.ishappy is zero or logical false, then POPULATE
%   was not able to obtain a happy result.
%
%   OP should be vectorized (i.e., accept a vector input), and ouput a vector of
%   the same length. Futhermore, OP may be a multi-valued function, in which
%   case it should accept a vector of length N and return a matrix of size NxM.
%
%   F.POPULATE(OP, VSCALE, HSCALE) enforces that the happiness check is relative
%   to the initial vertical scale VSCALE and horizontal scale HSCALE. These
%   values default to 0 and 1 respectively. During refinement, VSCALE updates
%   itself to be the largest magnitude values to which (each of the columns in)
%   OP evaluated to.
%
%   F.POPULATE(OP, VSCALE, HSCALE, PREF) enforces any additional preferences
%   specified in the preference structure PREF (see funcheb.PREF).
%
%   F.POPULATE(VALUES, ...) (or F.POPULATE({VALUES, COEFFS}, ...)) populates F
%   nonadaptively with the VALUES (and COEFFS) passed. These values are still
%   tested for happiness in the same way, but the length of the representation
%   is not reduced.
%
% See also FUNCHEB, FUNCHEB.pref, FUNCHEB.happinessCheck.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The (adaptive) FUNCHEB construction process is as follows:
%
%    --->[REFINE]      [values, flag] = pref.refinementStrategy(op, values, ...
%   |        |         pref). Allows refinements for: single sampling, 
%   |        |         resampling, and composition (see refine.m & compose.m).
%   |        v
%   |  [update VSCALE] VSCALE should only be computed from _sampled_ values, 
%   |        |         not extrapolated ones.
%   |        v
%   |   [EXTRAPOLATE]  Remove NaNs/Infs and (optionally) extrapolate endpoints.
%   |        |
%   |        v
%   | [compute COEFFS] COEFFS = chebpoly(VALUES)
%   |        |
%   |        v
%    -<--[ISHAPPY?]    [ISHAPPY, EPSLEVEL, CUTOFF] = pref.happinessCheck(f, op,
%            |         pref). Default check calls classicCheck() (previously 
%            | yes     'simplify.m') and sampleTest().
%            v
%      [alias COEFFS]  COEFFS = alias(COEFFS, TAIL)
%            |
%            v
%     [compute VALUES] VALUES = chebpolyval(COEFFS)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse inputs:
if ( nargin < 3 || isempty(vscale) )
    vscale = 0;
end
if ( nargin < 4 || isempty(hscale) )
    f.hscale = 1;
else
    f.hscale = hscale;
end
if ( nargin < 5 )
    pref = f.pref();
end

% Non adaptive construction. Values (and possibly coeffs) have been given.
if ( isnumeric(op) || iscell(op) )
    if ( isnumeric(op) )
        % OP is just the values.
        f.values = op;
        f.coeffs = f.chebpoly(op);
    else                 
        % OP is a cell {values, coeffs}
        f.values = op{1};
        f.coeffs = op{2};
    end
    % Update vscale:
    f.vscale = max(abs(f.values), [], 1);
    % Check for happiness: (no OP to compare against)
    [f.ishappy, f.epslevel] = happinessCheck(f, [], pref);
    return
end

% Initialise empty values to pass to refine:
values = [];

% Loop until ISHAPPY or GIVEUP:
while ( 1 )

    % Call the appropriate refinement routine:   
    [values, giveUp] = f.refine(op, values, pref);

    % We're giving up! :(
    if ( giveUp ) 
        break
    end    
    
    % Update vertical scale: (Only include sampled finite values)
    valuesTemp = values;
    valuesTemp(~isfinite(values)) = 0;
    vscale = max(vscale, max(abs(valuesTemp), [], 1));
    
    % Extrapolate out NaNs:
    [values, maskNaN, maskInf] = f.extrapolate(values);

    % Compute the Chebyshev coefficients:
    coeffs = f.chebpoly(values);
    
    % Check for happiness:
    f.values = values;
    f.coeffs = coeffs;
    f.vscale = vscale;
    [ishappy, epslevel, cutoff] = happinessCheck(f, op, pref); 
        
    % We're happy! :)
    if ( ishappy ) 
        coeffs = f.alias(coeffs, cutoff); % Alias the discarded coeffs.
        values = f.chebpolyval(coeffs);   % Compute values on this grid.
        break
    end
    
    % Replace any NaNs or Infs we may have extrapolated:
    values(maskNaN,:) = NaN;
    values(maskInf,:) = Inf;

end

% Update vertical scale one last time:
vscale = max(vscale, max(abs(values), [], 1));

% Assign to FUNCHEB object:
f.values = values;
f.coeffs = coeffs;
f.vscale = vscale;
f.ishappy = ishappy;
f.epslevel = epslevel;

% Extrapolate should have already dealt with NaNs and Infs if we were happy.
if ( ishappy )
    return
end

% Check for Infs: (if not happy)
if ( any(isinf(vscale)) )                       
    error('CHEBFUN:FUNCHEB:constructor:inf_blowup', ...
     'Function returned Inf when evaluated.')
end

end

