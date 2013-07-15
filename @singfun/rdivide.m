function s = rdivide(f, g)
%./   Divide SINGFUNS with SINGFUNS
%
%   This method will be called only if both F and G are SINGFUNS or at the 
%   most one of F and G is a scalar double.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%
% Check if inputs are other than SINGFUNS or doubles
if ( ~isa(f, 'singfun') && ~isa(f, 'double') )
    error( 'SINGFUN:rdivide:Input can only be a singfun or a double' )
end

if ( ~isa(g, 'singfun') && ~isa(g, 'double') )
    error( 'SINGFUN:rdivide:Input can only be singfun or a double' )
end

%%
% Reciprocal of a SINGFUN scaled by the double F.
if ( isa(f,'double') )    
    % convert f to a SINGFUN and call rdivide again.
    
    % Make a zero SINGFUN
    temp = singfun.zeroSingFun();        
    % Assign f to it's smooth part
    hscale = [];
    vscale  = [];
    smoothPrefs = chebtech.pref('tech', 'cheb1', 'extrapolate', false);
    temp.smoothPart = chebtech.constructor(f, vscale, hscale, smoothPrefs);
    % change f to a SINGFUN
    f = temp;
    % call SINGFUN.RDIVIDE again
    s = f./g;
    return
end

% Easy case of scalar division by a double
if ( isa(g,'double') )
    % copy the other input (a SINGUN) in the output
    s = f;
    % Divide the smooth part with the double and return
    s.smoothPart = (1/g) * f.smoothPart;
    return
end

%% SINGFUN./SINGFUN
% Division of two SINGFUNS
exponents = f.exponents - g.exponents;
if ( all(exponents > 0 ) ) 
    % division results in a smooth function with no singularities.
    s = singfun.zeroSingFun();
    smoothOp = @(x) feval(f, x)./feval(g, x);
    hscale = [];
    vscale  = [];
    smoothPrefs = chebtech.pref('tech', 'cheb1', 'extrapolate', false);
    s.smoothPart = chebtech.constructor(smoothOp, vscale, hscale, smoothPrefs);
else
    % Note: Exponents can be zero to generate a singular function. 
    % Example: f = 1; g = cos(pi/2*x) with trivial exponents. Then
    % s = f./g is singular with non trivial exponents. So the result of 
    % f./g in general is a generic SINGFUN with possibly non-trivial 
    % exponents    
    
    % factor out the known exponents
    op = singfun.singOp2SmoothOp( @(x) feval(f, x)./feval(g, x), exponents );
    % construct the singfun
    s = singfun( @(x) op(x), [], {'sing', 'sing'}, singfun.pref);    
    % add exponents
    s.exponents = s.exponents + exponents;
end

end