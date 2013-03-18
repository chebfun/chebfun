function [values, giveUp] = refine(op, values, pref)
%REFINE Refinement method for CHEBTECH1 construction.
%   VALUES = REFINE(OP, VALUES, PREF) determines the new VALUES of the operator
%   OP to be checked for happiness in the CHEBTECH1 construction process. The
%   exact procedure used is determined by PREF.CHEBTECH.REFINEMENTFUNCTION.
%
%   [VALUES, GIVEUP] = REFINE(OP, VALUES, PREF) returns also a binary GIVEUP
%   flag where TRUE means the refinement procedure has failed (typically when
%   the maximum number of samples, PREF.CHEBTECH.MAXSAMPLES, has been reached).
%
%   As opposed to CHEBTECH2, the only built-in refinement strategies is 'default'.
%   It makes use of grids of the form 2^(3:6 6:.5:16)+1 and resamples all of the 
%   values each time N is increased.
%
%   Alternative refinement strategies can be used by pasing a function handle in
%   the pref.chebtech.refinementStrategy field. The function handle should point
%   to a function with the template [VALUES, GIVEUP] = @(OP, VALUES, PREF) ...
%   which accepts a function handle OP, previously sampled values VALUES of OP
%   at a 1st-kind Chebyshev grid, and PREF, a preference structure containing
%   CHEBTECH1 preferences. It should return either a new set of VALUES (typically
%   on a finer grid) or set the GIVEUP flag to TRUE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Obtain some preferences:
if ( nargin < 3 )
    pref = chebtech.pref;
end
% No values were given:
if ( nargin < 2 )
    values = [];
end

% Grab the refinement function from the preferences:
refFunc = pref.chebtech.refinementFunction;

% Decide which refinement to use:
if ( ischar(refFunc) )
    % Re-sampling:
    [values, giveUp] = refineResampling(op, values, pref);
else
    % User defined refinement function:
    [values, giveUp] = refFunc(op, values, pref);
end
    
end

function [values, giveUp] = refineResampling(op, values, pref)
%REFINERESAMPLING  Default refinement function for resampling scheme.

    if ( isempty(values) )
        % Choose initial n based upon minSamples:
        n = 2^ceil(log2(pref.chebtech.minSamples) - 1) + 1;
    else
        % (Approximately) powers of sqrt(2):
        pow = log2(size(values, 1) - 1);
        if ( pow == floor(pow) && pow > 5 )
            n = round(2^(floor(pow) + .5)) + 1;
            n = n - mod(n, 2) + 1;
        else
            n = 2^(floor(pow)+1) + 1;
        end
    end
    
    % n is too large!
    if ( n > pref.chebtech.maxSamples )
        giveUp = true;
        return
    else
        giveUp = false;
    end
    
    % 1st-kind Chebyshev grid
    x = chebtech1.chebpts(n);

    values = feval(op, x);

end
