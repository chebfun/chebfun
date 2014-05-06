function val = feval(f, x)
%FEVAL   Evaluate a DELTAFUN.
%   FEVAL(F, X) evaluates the DELTAFUN F at the given points X. If the point of
%   evaluation has a non-trivial delta function, NaN is returned.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Trivial cases:
if ( isempty(x) )
    val = x;
    return
end
if ( isempty(f) )
    val = zeros(size(x));
    return
end

% Evaluate the smooth part.
val = feval(f.funPart, x);

% Point values of distributions do not make sense mathematically, so return NaN.
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
            idx = abs(x - deltaLoc(i))./deltaLoc(i) < proximityTol;
        end
        val(idx) = NaN;
    end
end
