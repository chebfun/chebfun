function [values, coeffs, vscale, ishappy, epslevel] = ...
    constructor(op, vscale, pref)
%CONSTRUCTOR    Construction process for FUNCHEB2 class.
%   CONSTRUCTOR(OP) returns the values of the function OP evaluated on a
%   sufficiently fine Chebyshev grid (of 2nd-kind points) such that the
%   Chebyshev interpolant through this data is deemed 'happy' (see
%   HAPPINESSCHECK.m). Essentially this means that such an interpolant is a
%   sufficiently accurate approximation to OP. OP should be vectorized (i.e.,
%   accept a vector input), and ouput a vector of the same length
%
%   OP may be a multi-valued function, in which case it should accept a vector
%   of length N and return a matrix of size NxM.
%
%   CONSTRUCTOR(OP, VSCALE) enforces that the happiness check is relative to the
%   vertical scale VSCALE, and CONSTRUCTOR(OP, VSCALE, PREF) enforces any
%   additional preferences specified in the preference structure PREF (see
%   funcheb2.PREF).
%
%   [VALUES, COEFFS] = CONSTRUCTOR(OP) returns both the values VALUES of the
%   interpolant and the coefficients COEFFS of the 1st-kind Chebyshev expansion
%   on the polynomial interpolant through this data.
%
%   [VALUES, COEFFS, VSCALE, ISHAPPY, EPSLEVEL] = CONSTRUCTOR(OP) additionally
%   returns an updated VSCALE, which is the largest magnitude at which OP was
%   sampled, ISHAPPY, which denotes whether the FUNCHEB2 is happy, and EPSLEVEL
%   which is the tolerance to which the FUNCHEB2 was deemed happy. If ISHAPPY is
%   zero or logical false, then CONSTRUCTOR was not able to obtain a happy
%   result. See HAPPINESSCHECK.m for further documentation on 'happiness'.
%
% See also FUNCHEB2, FUNCHEB2.pref, FUNCHEB2.happinessCheck.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The FUNCHEB2 construction process is as follows:
%
%    --->[REFINE]      [values, flag] = pref.funcheb2.refinementStrategy(op, ..
%   |        |         values, pref). Allows refinements for: single sampling, 
%   |        |         resampling, and composition (see refine.m & compose.m).
%   |        v
%   |  [update VSCALE] VSCALE should only be computed from _sampled_ values, 
%   |        |         not extrapolated ones.
%   |        v
%   |   [EXTRAPOLATE]  Remove NaNs/Infs and (optionally) extrapolate endpoints.
%   |        |
%   |        v
%   | [compute COEFFS] COEFFS = funcheb2.chebpoly(VALUES)
%   |        |
%   |        v
%    -<--[ISHAPPY?]    [ISHAPPY, EPSLEVEL, CUTOFF] = pref.happinessCheck(,
%     no     |         values coeffs, vscale, pref). Default check calls
%            | yes     funcheb2.classicCheck (prev 'simplify.m') and sampleTest.
%            v
%      [alias COEFFS]  COEFFS = funcheb2.alias(COEFFS, TAIL)
%            |
%            v
%     [compute VALUES] VALUES = funcheb2.chebpoly(COEFFS)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse inputs:
if ( nargin < 2 || isempty(vscale) )
    vscale = 0;
end
if ( nargin < 3 )
    pref = funcheb2.pref;
end

% Initialise empty values to pass to refine:
values = [];

% Loop until ISHAPPY or GIVEUP:
while ( 1 )

    % Call the appropriate refinement routine:   
    [values, giveUp] = funcheb2.refine(op, values, pref);

    % We're giving up! :(
    if ( giveUp ) 
        break
    end    
    
    % Update vertical scale: (Only include sampled finite values)
    valuesTemp = values;
    valuesTemp(~isfinite(values)) = 0;
    vscale = max(vscale, norm(valuesTemp(:), inf));
    
    % Extrapolate out NaNs:
    [values, maskNaN, maskInf] = funcheb2.extrapolate(values);

    % Compute the Chebyshev coefficients:
    coeffs = funcheb2.chebpoly(values);
    
    % Check for happiness:
    [ishappy, epslevel, cutoff] = ...
        funcheb2.happinessCheck(op, values, coeffs, vscale, pref); 
        
    % We're happy! :)
    if ( ishappy ) 
        coeffs = funcheb2.alias(coeffs, cutoff); % Alias the discarded coeffs.
        values = funcheb2.chebpolyval(coeffs);   % Compute values on this grid.
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
    error('CHEBFUN:FUNCHEB2:constructor:inf_blowup', ...
     'Function returned Inf when evaluated.')
end

% Check for NaNs: (if not happy)
if ( pref.funcheb2.extrapolate )
    % Check for NaNs in interior only:
    if ( any(any(isnan(values(2:end-1,:)))))
        error('CHEBFUN:FUNCHEB2:constructor:naneval', ...
            'Function returned NaN when evaluated.')
    end
    % We make sure to return something sensible (i.e., not a NaN) at +1 and -1.
    valuesTemp([1, end], :) = NaN;
    valuesTemp = funcheb2.extrapolate(valuesTemp);
    values([1, end],:) = valuesTemp([1, end],:);
else
    % Here we throw an error if NaNs were encountered anywhere.
    if ( any(isnan(values(:))) )
        error('CHEBFUN:FUNCHEB2:constructor:naneval2', ...
            'Function returned NaN when evaluated.')
    end
end

end

