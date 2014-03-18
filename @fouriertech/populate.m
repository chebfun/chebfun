function f = populate(f, op, vscale, hscale, pref)

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( (nargin < 3) || isempty(vscale) )
    vscale = 0;
end
if ( (nargin < 4) || isempty(hscale) )
    f.hscale = 1;
else
    f.hscale = hscale;
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Non-adaptive construction. %%%%%%%%%%%%%%%%%%%%%%%%%%
% Values (and possibly coefficients) have been given.
if ( isnumeric(op) || iscell(op) )
    if ( isnumeric(op) )
        % OP is just the values.
        f.values = op;
        f.isReal = isreal(op);
        f.coeffs = f.vals2coeffs(op);
    else                 
        % OP is a cell {values, coeffs}
        f.values = op{1};
        f.isReal = isreal(op);
        f.coeffs = op{2};
        if ( isempty(f.values) )
            f.values = f.coeffs2vals(f.coeffs);
            if max(abs(imag(f.values))) > 1e3*eps
                f.isReal = false;
            else
                f.isReal = true;
            end
        end
    end
    
    % Update vscale:
    f.vscale = max(abs(f.values), [], 1);
    
    % We're always happy if given discrete data:
    f.ishappy = true;
    
    % TODO: Is this the correct vscale?
    f.epslevel = max(eps(max(f.vscale)) + 0*f.vscale, eps);

    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Adaptive construction. %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise empty values to pass to refine:
f.values = [];

% Check a random value of op in (-pi,pi) to see if the result is complex.
f.isReal = isreal(feval(op,pi*(2*rand-1)));

% Loop until ISHAPPY or GIVEUP:
while ( 1 )

    % Call the appropriate refinement routine: (in PREF.REFINEMENTFUNCTION)
    [f.values, giveUp] = f.refine(op, f.values, pref);

    % We're giving up! :(
    if ( giveUp ) 
        break
    end    
    
    % Update vertical scale:
    f.vscale = max(vscale, max(abs(f.values), [], 1));

    
    % Compute the Fourier coefficients:
    f.coeffs = f.vals2coeffs(f.values);
    
    % Check for happiness:
    [ishappy, epslevel, cutoff] = happinessCheck(f, op, pref); 
        
    % We're happy! :)
    if ( ishappy ) 
        f.coeffs = f.alias(f.coeffs, cutoff); % Alias the discarded coefficients.
        f.values = f.coeffs2vals(f.coeffs);   % Compute values on this grid.
        break
    end

end

% This may not be the best place to do this.
if f.isReal
    f.values = real(f.values);
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Assign to FOURIERTECH object. %%%%%%%%%%%%%%%%%%%%%%%%%%
f.ishappy = ishappy;
f.epslevel = epslevel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ouput. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( ishappy )
    % We're done, and can return.
    return
end

end
