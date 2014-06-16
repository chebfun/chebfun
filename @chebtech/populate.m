function [f, values] = populate(f, op, vscale, hscale, pref)
%POPULATE   Populate a CHEBTECH class with values.
%   F = F.POPULATE(OP) returns a CHEBTECH representation populated with values
%   VALUES of the function OP evaluated on a Chebyshev grid. The fields
%   F.ISHAPPY and F.EPSLEVEL indicate whether the representation is deemed
%   'happy' and to what accuracy (see HAPPINESSCHECK.m). Essentially this means
%   that such an interpolant is a sufficiently accurate (i.e., to a relative
%   accuracy of F.EPSLEVEL) approximation to OP. If F.ISHAPPY is FALSE, then
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

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
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
%    -<--[ISHAPPY?]    [ISHAPPY, EPSLEVEL, CUTOFF] = PREF.HAPPINESSCHECK(F, OP,
%     no     |         PREF). Default calls CLASSICCHECK() and SAMPLETEST().
%            | yes     
%            v
%      [alias COEFFS]  COEFFS = ALIAS(COEFFS, CUTOFF)
%            |
%            v
%     [compute VALUES] VALUES = COEFFS2VALS(COEFFS)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse inputs:
if ( (nargin < 3) || isempty(vscale) )
    vscale = 0;
end
if ( (nargin < 4) || isempty(hscale) )
    f.hscale = 1;
else
    f.hscale = hscale;
end
if ( nargin < 5 )
    pref = chebtech.techPref();
end

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
    
    % Update vscale:
    f.vscale = getvscl(f);
    
    % We're always happy if given discrete data:
    f.ishappy = true;
    
    % Scale the epslevel relative to the largest column:
    vscale = f.vscale;
    f.epslevel = 10*eps(max(f.vscale));
    vscale(vscale <= f.epslevel) = 1;
    f.epslevel = f.epslevel./vscale;

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
    vscale = max(vscale, max(abs(valuesTemp)));
    
    % Extrapolate out NaNs:
    [values, maskNaN, maskInf] = extrapolate(f, values);

    % Compute the Chebyshev coefficients:
    coeffs = f.vals2coeffs(values);
    
    % Check for happiness:
    f.coeffs = coeffs;
    f.vscale = vscale;
    [ishappy, epslevel, cutoff] = happinessCheck(f, op, values, pref); 
        
    if ( ishappy ) % We're happy! :)
        % Alias the discarded coefficients:
        coeffs = f.alias(coeffs, cutoff);  
        break
    end
    
    % Replace any NaNs or Infs we may have extrapolated:
    values(maskNaN,:) = NaN;
    values(maskInf,:) = Inf;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update the vscale. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the 'true' vscale (as defined in CHEBTECH classdef):
vscaleOut = max(abs(values), [], 1);
% Update vertical scale one last time:
vscaleGlobal = max(vscale, vscaleOut);

% Output the 'true' vscale (i.e., the max of the stored values):
vscale = vscaleOut;

% Adjust the epslevel appropriately:
ind = vscaleOut < epslevel;
vscaleOut(ind) = epslevel(ind);
ind = vscaleGlobal < epslevel;
vscaleGlobal(ind) = epslevel(ind);
epslevel = epslevel.*vscaleGlobal./vscaleOut;

%%%%%%%%%%%%%%%%%%%%%%%%%% Assign to CHEBTECH object. %%%%%%%%%%%%%%%%%%%%%%%%%%
f.coeffs = coeffs;
f.vscale = vscale;
f.ishappy = ishappy;
f.epslevel = epslevel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ouput. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( ishappy )
    % We're done, and can return.
    f = simplify(f, f.epslevel/100);
    return
end

% [TODO]: Further decision needs to be made for whether the following lines
% should be deleted permanently. Currently they are commented out by KX.

% Check for Infs: (if not happy)
% if ( any(isinf(vscale)) )                       
%     error('CHEBFUN:CHEBTECH:populate:infBlowup', ...
%      'Function returned Inf when evaluated.')
% end

end
