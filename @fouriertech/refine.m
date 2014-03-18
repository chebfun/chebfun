function [values, giveUp] = refine(op, values, pref)
%REFINE   Refinement method for CHEBTECH2 construction.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Obtain some preferences:
% if ( nargin < 3 )
    pref = chebtech.techPref();
% end

% No values were given:
if ( nargin < 2 )
    values = [];
end

% Grab the refinement function from the preferences:
% refFunc = pref.refinementFunction;
% TODO: Implement other options.
refFunc = 'resampling';

% Decide which refinement to use:
if ( strcmpi(refFunc, 'nested') )
    % Nested ('single') sampling:
    [values, giveUp] = refineNested(op, values, pref);
elseif ( strcmpi(refFunc, 'resampling') )
    % Double sampling:
    [values, giveUp] = refineResampling(op, values, pref);
else
    % User defined refinement function:
%     [values, giveUp] = refFunc(op, values, pref);
end
    
end

function [values, giveUp] = refineResampling(op, values, pref)
%REFINERESAMPLING   Default refinement function for resampling scheme.

    if ( isempty(values) )
        % Choose initial n based upon minPoints:
        n = 2^ceil(log2(pref.minPoints - 1)) + 1;
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
    if ( n > pref.maxPoints )
        giveUp = true;
        return
    else
        giveUp = false;
    end
   
    % 2nd-kind Chebyshev grid:
%     x = chebtech2.chebpts(n);
    x = fourierpts(n);

    % Evaluate the operator:
    if ( pref.extrapolate )
        valuesTemp = feval(op, x(2:n-1));
        nans = NaN(1, size(valuesTemp, 2));
        values = [ nans; valuesTemp; nans ];
    else
        values = feval(op, x);
    end

end

function [values, giveUp] = refineNested(op, values, pref)
%REFINENESTED  Default refinement function for single ('nested') sampling.

    if ( isempty(values) )
        % The first time we are called, there are no values
        % and REFINENESTED is the same as REFINERESAMPLING.
        [values, giveUp] = refineResampling(op, values, pref);

    else
    
        % Compute new n by doubling (we must do this when not resampling).
        n = 2*size(values, 1) - 1;
        
        % n is too large:
        if ( n > pref.maxPoints )
            giveUp = true;
            return
        else
            giveUp = false;
        end
        
        % 2nd-kind Chebyshev grid:
        x = chebtech2.chebpts(n);
        % Take every 2nd entry:
        x = x(2:2:end-1);

        % Shift the stored values:
        values(1:2:n,:) = values;
        % Compute and insert new ones:
        values(2:2:end-1,:) = feval(op, x);

    end
end
