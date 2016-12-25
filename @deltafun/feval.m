function val = feval(f, x, varargin)
%FEVAL   Evaluate a DELTAFUN.
%   FEVAL(F, X) evaluates the DELTAFUN F at the given points X. If the point of
%   evaluation has a non-trivial delta function but no higher order delta 
%   functions, a signed infinity is returned. If howerver, there are higher
%   order delta functions present as well, a NaN is returned.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.


%% Trivial cases:
if ( isempty(f) )
    val = [];
    return
elseif ( isempty(x) )
    % Return empty matrix with dimensions of the appropriate size.
    val = zeros(size(x));
    return
end

%% LEFT AND RIGHT LIMITS: feval(f, x, 'left') and feval(f, x, 'right').
% Since delta functions are localised, in these limit cases, neglect delta 
% functions and call the feval() of the funPart:
if ( nargin > 2 && ~isempty(varargin) )
    val = feval(f.funPart, x, varargin);
    return
end

%% 
% Evaluate the smooth part.
val = feval(f.funPart, x);

pref = chebfunpref();
proximityTol = pref.deltaPrefs.proximityTol;        

% Make sure there are no trivial delta functions:
f = simplifyDeltas(f);
if ( isa(f, 'deltafun') )
    deltaLoc = f.deltaLoc;
    for i = 1:length(deltaLoc)
        if ( deltaLoc(i) == 0 )
            % Avoid divide by zero:
            idx = abs(x - deltaLoc(i)) < proximityTol;
        else
            % Check the relative distance from delta function locations:
            idx = abs((x - deltaLoc(i))./deltaLoc(i)) < proximityTol;
        end
        
        if ( abs(f.deltaMag(1, i)) > 0 )
            % If there is a delta function, assign a signed infinity:
            val(idx) = Inf*sign(f.deltaMag(1, i));
        else
            % If there is no delta function but higher order, assign NaNs:
            val(idx) = NaN;
        end
    end
end
end
