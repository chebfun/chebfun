function h = times(f, g)
%.*   Multiply SINGFUNS with SINGFUNS and SMOOTHFUNS
%   F.*G multiplies SINGFUN objects F and G or a SINGFUN by a scalar/SMOOTHFUN
%   if either F or G is a scalar/SMOOTHFUN.
%
% See also LDIVIDE, RDIVIDE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%
% Empty arguments:
if ( isempty(f) || isempty(g) )
    h = singfun();
    return
end

%%
% Check if inputs are other than SINGFUNS, SMOOTHFUNS or doubles:
if ( (~isa(f, 'singfun') && ~isa(f, 'smoothfun') && ~isa(f, 'double')) || ...
     (~isa(g, 'singfun') && ~isa(g, 'smoothfun') && ~isa(g, 'double')) )
    error('CHEBFUN:SINGFUN:times:input', ...
        'Input can only be a singfun, a smoothfun or a double')
end

% Multiplication by a scalar or a SMOOTHFUN:
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

h = cancelExponents(h);

%% Simplify and replace the boundary roots:
h = simplify(h);

%%
% Check if after multiplication h has become smooth:
if ( issmooth(h) )
    h = h.smoothPart;
end

end
