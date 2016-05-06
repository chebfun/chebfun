function [values, giveUp] = refine(op, values, pref)
%REFINE   Refinement method for CHEBTECH1 construction.
%   VALUES = REFINE(OP, VALUES, PREF) determines the new VALUES of the operator
%   OP to be checked for happiness in the CHEBTECH1 construction process. The
%   exact procedure used is determined by PREF.REFINEMENTFUNCTION.
%
%   [VALUES, GIVEUP] = REFINE(OP, VALUES, PREF) returns also a binary GIVEUP
%   flag where TRUE means the refinement procedure has failed (typically when
%   the maximum number of points, PREF.MAXLENGTH, has been reached).
%
%   The two built-in refinement strategies are 'NESTED' and 'RESAMPLING'. The
%   former makes use of the nested property of the 1st-kind grid by taking N
%   (the number of points) to be 17*3^(0:7) and doesn't resample previously
%   evaluated values. The latter uses grids of the form 2^(4:6 6:.5:16) + 1 and
%   resamples all of the values each time N is increased. The 'RESAMPLING'
%   option should be used for functions which are not sampleable, for example,
%   anything that depends on the length of the input to OP.
%
%   Alternative refinement strategies can be used if the
%   PREF.REFINEMENTFUNCTION field is a function handle. The function handle
%   should point to a function with the template [VALUES, GIVEUP] = @(OP,
%   VALUES, PREF) ... which accepts a function handle OP, previously sampled
%   values VALUES of OP at a 1st-kind Chebyshev grid, and PREF, a preference
%   structure containing CHEBTECH preferences. It should return either a new
%   set of VALUES (typically on a finer grid) or set the GIVEUP flag to TRUE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Obtain some preferences:
if ( nargin < 3 )
    pref = chebtech.techPref();
end

% No values were given:
if ( nargin < 2 )
    values = [];
end

% Grab the refinement function from the preferences:
refFunc = pref.refinementFunction;

% Decide which refinement to use:
if ( strcmpi(refFunc, 'nested') )
    % Nested ('single') sampling:
    [values, giveUp] = refineNested(op, values, pref);
elseif ( strcmpi(refFunc, 'resampling') )
    % Double sampling:
    [values, giveUp] = refineResampling(op, values, pref);
else
    % User defined refinement function:
    [values, giveUp] = refFunc(op, values, pref);
end

% Ensure that doubles are returned:
values = double(values);
    
end

function [values, giveUp] = refineResampling(op, values, pref)
%REFINERESAMPLING   Default refinement function for resampling scheme.

    if ( isempty(values) )
        % Choose initial n based upon minSamples:
        n = 2^ceil(log2(pref.minSamples - 1)) + 1;
    else
        % (Approximately) powers of sqrt(2):
        pow = log2(size(values, 1) - 1);
        if ( (pow == floor(pow)) && (pow > 5) )
            n = round(2^(floor(pow) + .5)) + 1;
            n = n - mod(n, 2) + 1;
        else
            n = 2^(floor(pow) + 1) + 1;
        end
    end
    
    % n is too large:
    if ( n > pref.maxLength )
        % Don't give up if we haven't sampled at least once.
        if ( isempty(values) )
            n = pref.maxLength;
            giveUp = false;
        else
            giveUp = true;
            return
        end
    else
        giveUp = false;
    end
    
    % 1st-kind Chebyshev grid:
    x = chebtech1.chebpts(n);

    values = feval(op, x);

end

function [values, giveUp] = refineNested(op, values, pref)
%REFINENESTED  Default refinement function for single ('nested') sampling.

    if ( isempty(values) )
        % The first time we are called, there are no values
        % and REFINENESTED is the same as REFINERESAMPLING.
        [values, giveUp] = refineResampling(op, values, pref);

    else
    
        % Grab the number of the points used in last iteration:
        n = size(values, 1);
        
        if ( ( n < pref.maxLength ) && ( 3*n > pref.maxLength ) )
            
            % Set N as maxLength to make best of the largest number of points 
            % allowed:
            n = pref.maxLength;
            
            % 1st-kind Chebyshev grid:
            x = chebtech1.chebpts(n);
            
            % Evaluate OP at X:
            values = feval(op, x);
            giveUp = false;

        elseif ( n < pref.maxLength )
            
            % triple the points:
            n = 3*n;
            
            % 1st-kind Chebyshev grid:
            x = chebtech1.chebpts(n);
            
            % Re-group the points:
            x1 = x(1:3:end-2);
            x3 = x(3:3:end);
            
            % Copy of the sampled function values:
            oldValues = values;
            
            % Compute and insert new ones:
            values(1:3:n,:) = feval(op, x1);
            values(3:3:n,:) = feval(op, x3);
            
            % Re-distribute the stored values:
            values(2:3:n,:) = oldValues;
        
            giveUp = false;
        else
            
            % For any other cases, give up:
            giveUp = true;
            return
        end
        
    end
end
