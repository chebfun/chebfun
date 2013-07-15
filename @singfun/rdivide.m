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
s = singfun;
if ( all(f.exponents - g.exponents >= 0 ) )
    % division results in a smooth function with no singularity
    s = singfun( @(x) feval(./feval(g.smoothPart, x), [], {'sing', 'sing'}, singfun.pref);      
% multiply the smooth parts
s.smoothPart = (f.smoothPart)./(g.smoothPart);
% add the exponents
s.exponents = f.exponents - g.exponents;

%%
% Check if after division the type of singularity has changed or if 
% it can be removed.
% [TODO]: Since exponents are negative, it's impossible to remove a 
% singularity after mutiplying two SINGFUNS?
tol = singfun.pref.singfun.eps;
% loop through each end
for k = 1:2
    if ( s.exponents(k) < -100*tol )
        s.isSingEnd(k) = 1;
        if ( abs(s.exponents(k) - round(s.exponents(k))) < 100*tol )
            s.singType{k} = 'pole';
        else
            s.singType{k} = 'sing';
        end
    else
        s.isSingEnd(k) = 0;
        s.singType{k} = 'none';
    end
end

end