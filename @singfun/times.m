function h = times(f, g)
%.*   Multiply SINGFUNS with SINGFUNS
%   F.*G multiplies SINGFUN objects F and G or a SINGFUN by a scalar if either
%   F or G is a scalar.

% See also LDIVIDE, RDIVIDE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note: This method will be called only if both F and G are SINGFUNS or at the
% most one of F and G is a scalar double.

% NH: Should also support SINGFUN.*SMOOTHFUN.

%% Trivial cases:

% Check if inputs are other than SINGFUNS or doubles:
if ( (~isa(f, 'singfun') && ~isa(f, 'double')) || ...
        (~isa(g, 'singfun') && ~isa(g, 'double')) )
    error( 'SINGFUN:times:Input can only be a singfun or a double' )
end

% Scalar multiplication cases
if ( isa(f, 'double') )
    % Copy the other input (a SINGFUN) in the output:
    h = g;
    % Multiply the smooth part with the double and return:
    h.smoothPart = f * g.smoothPart;
    return
elseif ( isa(g, 'double') )
    % Copy the other input (a SINGFUN) in the output:
    h = f;
    % Multiply the smooth part with the double and return
    h.smoothPart = g * f.smoothPart;
    return
end

%% Multiplication of two SINGFUNS:

% Initialise the output SINGFUN:
h = singfun();
% Multiply the smooth parts:
h.smoothPart = (f.smoothPart) .* (g.smoothPart);
% Add the exponents:
h.exponents = f.exponents + g.exponents;

%%
% Check if after multiplication the singularity can be removed.
tol = singfun.pref.singfun.exponentTol;
h.exponents(abs(h.exponents) < tol) = 0;
% NH: h = h.smoothPart; ?

end
