function f = fracCalc(f, m)
%FRACCALC   Differentiation or integral of a BNDFUN for a fractional order.
%  FRACCALC(F, M) returns the fractional derivative or the integral of F 
%  which is defined on a finite domain.
%
% See also CUMSUM, DIFF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Trivial case:
if ( isempty(f) )
    return
end

% Get the exponents of u
exps = get(f.onefun ,'exponents');


    % integrand of the operator
    

newexps = exps;
newexps(1) = exps(1)+alpha;
newends = ends;



end

function y = integrand(x)

% fractional kernel
k = @(x,s) (x-s).^(alpha-1);

if any(x == ends(1))
    y = chebfun(0,[ends(1),ends(1)]);
elseif any(x == ends(2:end))
    y = chebfun(NaN,[ends(1),x]);
else
    % y = chebfun(@(s) feval(u,s).*k(x,s),[ends(ends<x) x],'exps',[exps(1) alpha-1],'scale',u.scl)
    
    % playing with piecewise chebfuns
    newends = [ends(ends<x) x];
    tmpexps = [];
    for l = 1:length(newends)-1
        tmpexps = [tmpexps exps(l,1) 0];
    end
    tmpexps(end) = alpha-1;
    y = chebfun(@(s) feval(u,s).*k(x,s),newends,'exps',tmpexps,'scale',u.scl,'extrapolate',true);
    
end

end