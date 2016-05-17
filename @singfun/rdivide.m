function s = rdivide(f, g)
%./   Divide SINGFUNS with SINGFUNS
%   RDIVIDE(F, G) computes the pointwise division F./G. 
%   This method will be called only if both F and G are SINGFUNS or at the 
%   most one of F and G is a scalar double or a SMOOTHFUN. It is also 
%   assumed that G has no roots in the open interval (-1, 1). 
%
% See also LDIVIDE, TIMES.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%
% Case of empty arguments.
if ( isempty(f) || isempty(g) )
    % Return an empty SINGFUN:
    s = singfun();
    return
end

% Check if inputs are other than SINGFUNs, SMOOTHFUNs or doubles.
if ( (~isa(f, 'singfun') && ~isa(f, 'smoothfun') && ~isa(f, 'double')) || ...
     (~isa(g, 'singfun') && ~isa(g, 'smoothfun') && ~isa(g, 'double')) )   
    error('CHEBFUN:SINGFUN:rdivide:rdivide', ...
        'Input can only be a singfun, a smoothfun or a double')
end

%% SINGFUN./DOUBLE or SINGFUN./SMOOTHFUN
if ( isa(g,'double') || isa(g, 'smoothfun') )
    % Copy the other input (a SINGUN) in the output:
    s = f;
    % Divide the smooth part with the double g and return:
    s.smoothPart = (1./g) .* f.smoothPart;
    if ( issmooth(s) )
        s = s.smoothPart;
    end
    return
end

%% Reciprocal: DOUBLE./SINGFUN
if ( isa(f, 'double') )    
    % Make a SMOOTHFUN of the double f:
    f = g.smoothPart.make(f);
    % Convert f to a SINGFUN:
    f = singfun.smoothFun2SingFun(f);
end

%% SMOOTHFUN./SINGFUN
if ( isa(f, 'smoothfun') )    
    % Convert f to a SINGFUN and call rdivide again.
    f = singfun.smoothFun2SingFun(f);
end

%% SINGFUN./SINGFUN
% Note: Exponents of f and g can all be zero to generate a singular function.
% Example: f = 1; g = cos(pi/2*x) with trivial exponents. Then s = f./g is 
% singular with non trivial exponents. So the result of f./g in general is a
% generic SINGFUN with possibly non-trivial exponents.

% Note: Once we reach here, both f and g are SINGFUN objects.

% Extract boundary roots:
g = extractBoundaryRoots(g);

% Grab the boundary values of the smooth part of G:
boundaryValues = [get(g.smoothPart, 'lval') get(g.smoothPart, 'rval')];

% Set a tolerance:
tol = 1e2*get(g, 'vscale')*eps;

if ( all(abs(boundaryValues) > tol) )
    % No vanishing boundary values, then take advantage of the information
    % we know about the exponents:

    pref.extrapolate = 1;
    h = f.constructSmoothPart(@(x) feval(f.smoothPart, x)./ ...
        feval(g.smoothPart, x), [], pref);
    data.exponents = f.exponents - g.exponents;
    s = singfun(h, data);

else
    % Construct the SINGFUN by a direct call to the constructor:
    s = singfun(@(x) feval(f, x)./feval(g, x));
end

%% Simplify and replace the boundary roots:
s = cancelExponents(s);
s = simplify(s);

%% 
% Check if after division s has become smooth:
if ( issmooth(s) )
    s = s.smoothPart;
end

end
