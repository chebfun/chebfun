function s = times(f,g)
%.*   Multiply SINGFUNS with SINGFUNS

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note: This method will be called only if both F and G are SINGFUNS or at the
% most one of F and G is a scalar double.

%% Trivial cases:

% Check if inputs are other than SINGFUNS or doubles:
if ( ~isa(f, 'singfun') && ~isa(f, 'double') ) || ...
        ( ~isa(g, 'singfun') && ~isa(g, 'double') )
    error( 'SINGFUN:times:Input can only be a singfun or a double' )
end

% Scalar multiplication cases
if ( isa(f, 'double') )
    % Copy the other input (a SINGUN) in the output:
    s = g;
    % Multiply the smooth part with the double and return:
    s.smoothPart = f * g.smoothPart;
    return
elseif ( isa(g, 'double') )
    % Copy the other input (a SINGUN) in the output:
    s = f;
    % Multiply the smooth part with the double and return
    s.smoothPart = g * f.smoothPart;
    return
end

%% Mutliplication of two SINGFUNS:

% Initialise the output SINGFUN:
s = singfun;
% Multiply the smooth parts:
s.smoothPart = (f.smoothPart).*(g.smoothPart);
% Add the exponents:
s.exponents = f.exponents + g.exponents;

% Check if after multiplication the type of singularity has changed or if it can
% be removed.

% [TODO]: Since exponents are negative, it's impossible to remove a singularity
% after mutiplying two SINGFUNS?

tol = singfun.pref.singfun.eps;

% Loop through each end:
for k = 1:2
    if ( s.exponents(k) < -100*tol )
        if ( abs(s.exponents(k) - round(s.exponents(k))) < 100*tol )
            s.singType{k} = 'pole';
            s.exponents(k) = round(s.exponents(k));
        else
            s.singType{k} = 'sing';
        end
    else
        if ( s.exponents(k) < 100*tol )
            s.singType{k} = 'none';
            s.exponents(k) = 0;
        else
            s.singType{k} = 'root';
        end
    end
end

end
