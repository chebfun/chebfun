function s = rdivide(f,g)
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
tol = singfun.pref.singfun.eps;
% Reciprocal of a SINGFUN scaled by the double F.
if ( isa(f,'double') )    
    % take the scaled reciprocal of the SINGFUN G.
    if ( any(g.exponents < -100*tol) )
        % if g is singular the reciprocal of g will be a smooth function 
        % with no singularities
        s = singfun.zeroSingFun();
        smoothOp = @(x) 1./feval(g, x);
        hscale = [];
        vscale  = [];
        smoothPrefs = chebtech.pref('tech', 'cheb1', 'extrapolate', false);
        s.smoothPart = chebtech.constructor(smoothOp, vscale, hscale, smoothPrefs);
    else
        % g has trivial exponents, but the reciporcal of g might be
        % singular, for example, if the smooth part of g vanishes at an 
        % end point. So, call the SINGFUN constructor.
        s = singfun( @(x) 1./feval(g.smoothPart, x), [], {'sing', 'sing'}, singfun.pref);      
    end
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

%%
% Division of two SINGFUNS
if ( all(f.exponents - g.exponents > 0 ) )
    % division results in a smooth function with no singularity
    s = singfun.zeroSingFun();
    smoothOp = @(x) feval(f, x)./feval(g, x);
    hscale = [];
    vscale  = [];
    smoothPrefs = chebtech.pref('tech', 'cheb1', 'extrapolate', false);
    s.smoothPart = chebtech.constructor(smoothOp, vscale, hscale, smoothPrefs);
else
    % The result of f./g is a generic SINGFUN with possibly non-trivial exponents
    s = singfun( @(x) feval(f, x)./feval(g, x), [], {'sing', 'sing'}, singfun.pref);    
end

%%
end