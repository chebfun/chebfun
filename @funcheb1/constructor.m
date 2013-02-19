function [values, coeffs, vscale, ishappy, epslevel] = ...
    constructor(op, vscale, pref)
%CONSTRUCTOR    Construction process for FUNCHEB1 class.
%   CONSTRUCTOR(OP) returns the values of the function OP evaluated on a
%   sufficiently fine Chebyshev grid (of 1st-kind points) such that the
%   Chebyshev interpolant through this data is deemed 'happy' (see 
%   HAPPINESSCHECK.m). Essentially this means that such an interpolant is a
%   sufficiently accurate approximation to OP. OP should be vectorized (i.e.,
%   accept a vector input), and ouput a vector of the same length.
%
%   OP may be a multi-valued function, in which case it should accept a vector
%   of length N and return a matrix of size NxM.
%
%   CONSTRUCTOR(OP, VSCALE) enforces that the happiness check is relative to the
%   vertical scale VSCALE, and CONSTRUCTOR(OP, VSCALE, PREF) enforces any
%   additional preferences specified in the preference structure PREF (see
%   funcheb1.PREF).
%
%   [VALUES, COEFFS] = CONSTRUCTOR(OP) returns both the values VALUES of the
%   interpolant and the coefficients COEFFS of the 1st-kind Chebyshev expansion
%   on the polynomial interpolant through this data.
%
%   [VALUES, COEFFS, VSCALE, ISHAPPY, EPSLEVEL] = CONSTRUCTOR(OP) additionally
%   returns an updated VSCALE, which is the largest magnitude at which OP was
%   sampled, ISHAPPY, which denotes whether the FUNCHEB1 is happy, and EPSLEVEL
%   which is the tolerance to which the FUNCHEB1 was deemed happy. If ISHAPPY is
%   zero or logical false, then CONSTRUCTOR was not able to obtain a happy
%   result. See HAPPINESSCHECK.m for further documentation on 'happiness'.
%
% See also FUNCHEB1, FUNCHEB1.pref, FUNCHEB1.happinessCheck.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The FUNCHEB1 construction process is as follows:
%
%    --->[REFINE]      [values, flag] = pref.funcheb1.refinementStrategy(op, ..
%   |        |         values, pref). Allows refinements for: single sampling, 
%   |        |         resampling, and composition (see refine.m & compose.m).
%   |        v
%   |  [update VSCALE] VSCALE should only be computed from _sampled_ values, 
%   |        |         not extrapolated ones.
%   |        v
%   |   [EXTRAPOLATE]  Remove NaNs/Infs by extrapolation.
%   |        |
%   |        v
%   | [compute COEFFS] COEFFS = funcheb1.chebpoly(VALUES)
%   |        |
%   |        v
%    -<--[ISHAPPY?]    [ISHAPPY, EPSLEVEL, CUTOFF] = pref.happinessCheck(,
%     no     |         values coeffs, vscale, pref). Default check calls 
%            | yes     funcheb1.classicCheck (prev 'simplify.m') and sampleTest.
%            v
%      [alias COEFFS]  COEFFS = funcheb1.alias(COEFFS, TAIL)
%            |
%            v
%     [compute VALUES] VALUES = funcheb1.chebpoly(COEFFS)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse inputs:
if ( nargin < 2 || isempty(vscale) )
    vscale = 0;
end
if ( nargin < 3 )
    pref = funcheb1.pref;
end

% Initialise empty values to pass to refine:
values = [];

% Loop until ISHAPPY or GIVEUP:
while ( 1 )

    % Call the appropriate refinement routine:   
    [values, giveUp] = funcheb1.refine(op, values, pref);

    % We're giving up! :(
    if ( giveUp )
        break
    end    

    % Update vertical scale (Only include sampled finite values):
    valuesTemp = values;
    valuesTemp(~isfinite(values)) = 0;
    vscale = max(vscale, norm(valuesTemp(:), inf));

    % Extrapolate out NaNs:
    [values, maskNaN, maskInf] = funcheb1.extrapolate(values);

    % Compute the Chebyshev coefficients:
    coeffs = funcheb1.chebpoly(values);
    
    % Check for happiness:
    [ishappy, epslevel, cutoff] = ...
        funcheb1.happinessCheck(op, values, coeffs, vscale, pref); 
        
    % We're happy! :)
    if ( ishappy ) 
        coeffs = funcheb1.alias(coeffs, cutoff); % Alias the discarded coeffs.
        values = funcheb1.chebpolyval(coeffs); % Compute values on this grid.
        break
    end
    
    % Replace any NaNs or Infs we may have extrapolated:
    values(maskNaN,:) = NaN;
    values(maskInf,:) = Inf;

end

% Update vertical scale one last time:
vscale = max( vscale, norm(values(:), inf) );

% Extrapolate should have already dealt with NaNs and Infs if happy.
if ( ishappy )
    return
end

% Check for Infs: (if not happy)
if ( any(isinf(vscale)) )
    error('CHEBFUN:FUNCHEB1:constructor:inf_blowup', ...
     'Function returned Inf when evaluated.')
end

% Check for NaNs: (if not happy)
if ( any( any( isnan( values(:) ) ) ) )
    error('CHEBFUN:FUNCHEB1:constructor:naneval', ...
     'Function returned NaN when evaluated.')
end

