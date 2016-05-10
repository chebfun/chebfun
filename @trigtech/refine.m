function [values, giveUp] = refine(op, values, pref)
%REFINE   Refinement method for TRIGTECH construction.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Obtain some preferences:
if ( nargin < 3 )
    pref = trigtech.techPref();
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
    error('CHEBFUN:TRIGTECH:refine', ...
          'No user defined refinement options allowed')
end
    
end

function [values, giveUp] = refineResampling(op, values, pref)
%REFINERESAMPLING   Default refinement function for resampling scheme.

    if ( isempty(values) )
        % Choose initial n based upon minSamples:
        n = 2^ceil(log2(pref.minSamples - 1));
    else
        % (Approximately) powers of sqrt(2):
        pow = log2(size(values, 1));
        if ( (pow == floor(pow)) && (pow > 5) )
            n = 3*2^(pow-1);
        else
            n = 2^(floor(pow) + 1);
        end
    end
    % n is too large:
    if ( n > pref.maxLength )
        giveUp = true;
        return
    else
        giveUp = false;
    end
   
    % [TODO]: Allow "first-kind" trigonometric points.
    x = trigpts(n);    

    % Evaluate the operator:
    values = feval(op, x);
    
    % Force the value at -1 to be equal to the value at 1, thus
    % enforcing symmetry.
    valRightBoundary = feval(op,1);
    
    values(1,:) = 0.5*(values(1,:) + valRightBoundary);
end

function [values, giveUp] = refineNested(op, values, pref)
%REFINENESTED   Default refinement function for single ('nested') sampling.

    if ( isempty(values) )
        % The first time we are called, there are no values
        % and REFINENESTED is the same as REFINERESAMPLING.
        [values, giveUp] = refineResampling(op, values, pref);

    else
    
        % Compute new n by doubling (we must do this when not resampling).
        n = 2*size(values, 1);
        
        % n is too large:
        if ( n > pref.maxLength )
            giveUp = true;
            return
        else
            giveUp = false;
        end
        
        % "2nd-kind" trigonometric points:
        x = trigpts(n);
        % Take every 2nd entry:
        x = x(2:2:n);

        % Shift the stored values:
        values(1:2:n,:) = values;
        % Compute and insert new ones:
        values(2:2:n,:) = feval(op, x);

    end
    
end
