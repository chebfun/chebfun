function s = times(f,g)
%.*   Multiply DELTAFUNS with DELTAFUNS

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note: This method will be called only if both F and G are DELTAFUNS or at the
% most one of F and G is a scalar double.

%% Empty case:
s = deltafun;
if ( isempty(f) || isempty(g) )
    return
end

% Trivial cases:
% Check if inputs are other than DELTAFUNS or doubles:
if ( (~isa(f, 'deltafun') && ~isa(f, 'double')) || ...
     (~isa(g, 'deltafun') && ~isa(g, 'double')) )
    error( 'DELTAFUN:times', 'Input can only be a DELTAFUN or a double' )
end

%% Multiplication by a scalar or a function:
if ( isa(f, 'double') || isa(f, 'smoothfun') )
    % Copy the other input (a SINGFUN) in the output:
    h = g;
    % Multiply the smooth parts and return:
    h.smoothPart = f .* g.smoothPart;
elseif ( isa(g, 'double') || isa(g, 'smoothfun') )
    % Copy the other input (a SINGFUN) in the output:
    h = f;
    % Multiply the smooth parts and return:
    h.smoothPart = g .* f.smoothPart;
end

%% Multiplication of two SINGFUNS:
if ( isa(f, 'singfun') && isa(g, 'singfun') )
    % Initialise the output SINGFUN:
    h = singfun();
    % Multiply the smooth parts:
    h.smoothPart = (f.smoothPart) .* (g.smoothPart);
    % Add the exponents:
    h.exponents = f.exponents + g.exponents;
end

%%
% Check if after multiplication h has become smooth:
if ( issmooth(h) )
    h = h.smoothPart;
end

end
