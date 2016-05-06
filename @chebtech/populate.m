function [f, values] = populate(f, op, data, pref)
%POPULATE   Populate a CHEBTECH class with values.
%   F = F.POPULATE(OP) returns a CHEBTECH representation populated with values
%   VALUES of the function OP evaluated on a Chebyshev grid. The field
%   F.ISHAPPY indicates whether the representation is deemed 'happy' (see
%   HAPPINESSCHECK.m). Essentially this means that such an interpolant is a
%   sufficiently accurate approximation to OP.  If F.ISHAPPY is FALSE, then
%   POPULATE was not able to obtain a happy result.
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
%   specified in the preference structure PREF (see CHEBTECH.TECHPREF).
%
%   F.POPULATE(VALUES, ...) (or F.POPULATE({VALUES, COEFFS}, ...)) populates F
%   non-adaptively with the VALUES (and COEFFS) passed. These values are still
%   tested for happiness in the same way as described above, but the length of
%   the representation is not altered.
%
% See also CHEBTECH, TECHPREF, HAPPINESSCHECK.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The (adaptive) CHEBTECH construction process is as follows:
%
%    --->[REFINE]      [values, flag] = pref.refinementFunction(op, values, ...
%   |        |         pref). Allows refinements for: nested sampling, 
%   |        |         resampling, and composition (see REFINE.m & COMPOSE.m).
%   |        v
%   |  [update VSCALE] VSCALE should only be computed from _sampled_ values, 
%   |        |         not extrapolated ones.
%   |        v
%   |   [EXTRAPOLATE]  Remove NaNs/Infs and (optionally) extrapolate endpoints.
%   |        |
%   |        v
%   | [compute COEFFS] COEFFS = VALS2COEFFS(VALUES)
%   |        |
%   |        v
%    -<--[ISHAPPY?]    [ISHAPPY, CUTOFF] = PREF.HAPPINESSCHECK(F, OP, PREF).
%     no     |         Default calls STANDARDCHECK() and SAMPLETEST().
%            | yes     
%            v
%      [alias COEFFS]  COEFFS = ALIAS(COEFFS, CUTOFF)
%            |
%            v
%     [compute VALUES] VALUES = COEFFS2VALS(COEFFS)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Non-adaptive construction. %%%%%%%%%%%%%%%%%%%%%%%%%%
% Values (and possibly coefficients) have been given.
if ( isnumeric(op) || iscell(op) )
    values = op;
    if ( isnumeric(op) )
        % OP is just the values.
        if ( all(isnan(op)) )
            values = op;
        else
            values = extrapolate(f, values); 
        end
        f.coeffs = f.vals2coeffs(values);
    else                 
        % OP is a cell {values, coeffs}
        f.coeffs = op{2};
    end

    % We're always happy if given discrete data:
    f.ishappy = true;

    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Adaptive construction. %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise empty values to pass to refine:
values = [];

% Loop until ISHAPPY or GIVEUP:
while ( 1 )

    % Call the appropriate refinement routine: (in PREF.REFINEMENTFUNCTION)
    [values, giveUp] = f.refine(op, values, pref);

    % We're giving up! :(
    if ( giveUp ) 
        break
    end    
    
    % Update vertical scale: (Only include sampled finite values)
    valuesTemp = values;
    valuesTemp(~isfinite(values)) = 0;
    data.vscale = max(data.vscale, max(abs(valuesTemp)));
    
    % Extrapolate out NaNs:
    %
    % TODO:  This will extrapolate out however many NaNs / Infs there are, even
    % if that number is larger than the number of function values.  This is
    % probably not good.
    [values, maskNaN, maskInf] = extrapolate(f, values);

    % Compute the Chebyshev coefficients:
    coeffs = f.vals2coeffs(values);
    
    % Check for happiness:
    f.coeffs = coeffs;
    [ishappy, cutoff] = happinessCheck(f, op, values, data, pref);
        
    if ( ishappy ) % We're happy! :)
        % disard unwanted coefficients
        f = prolong(f,cutoff);  
        break
    end
    
    % Replace any NaNs or Infs we may have extrapolated:
    values(maskNaN,:) = NaN;
    values(maskInf,:) = Inf;

end

%%%%%%%%%%%%%%%%%%%%%%%%%% Assign to CHEBTECH object. %%%%%%%%%%%%%%%%%%%%%%%%%%
f.ishappy = ishappy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ouput. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( ishappy )
    % We're done, and can return.
    return
end

end
