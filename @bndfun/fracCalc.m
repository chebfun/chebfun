function f = fracCalc(funs, fracM)
%FRACCALC   Differentiation or integral of a BNDFUN for a fractional order.
%  FRACCALC(F, M) returns the fractional derivative or the integral of F 
%  which is defined on a finite domain.
%
% See also CUMSUM, DIFF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Check if F is built on SINGFUN, if not cast its ONEFUN to a SINGFUN:
numFuns = numel(funs);

for k = 1:numFuns
    if ( ~issing(funs{k}) )
        funs{k}.onefun = singfun(funs{k}.onefun);
    end
end

exps = funs{end}.onefun.exponents;
if ( numFuns == 1 )
    exps(1) = exps(1) + fracM;
end

f = funs{end};
f.onefun = (diff(funs{end}.domain)/2)^fracM*(1/gamma(fracM))*...
    singfun(@(x) op(funs, x, fracM), exps);

end

%%
function y = op(funs, x, fracM)

% Get the exponents of G:
g = funs{end}.onefun;
exps = g.exponents;
sp = g.smoothPart;

oldSize = size(x);
x = x(:);
l = length(x);
y = zeros(l, 1);

numFuns = numel(funs);

maps = cell(1, numFuns);
for j = 1:numFuns
    maps{j} = funs{j}.mapping;
end

for k = 1:l
        
    % Compute the contribution from the previous pieces:
    I = 0;
    
    for j = 1:(numFuns-1)
        
        t = maps{j}.inv(maps{end}.for(x(k)));
        kernel = singfun(@(s) (t-s).^(fracM-1), [0 0]);
        Ij = sum(kernel.*funs{j}.onefun);
        
        I = I + Ij;
    end
    
    if ( x(k) == -1 )
        % When sampled at -1, the definite integral given above, i.e.
        % sum(integrand(f, x)) should be zero, which is equivalent to a zero
        % SINGFUN:
        y(k) = NaN;
    elseif ( x(k) == 1 )
        % Call SINGFUN constructor:
        h = singfun(sp, [exps(1) exps(2)+fracM-1]);
        y(k) = I + sum(h);
    else
        % Call SINGFUN constructor:
        rsp = restrict(sp, [-1 x(k)]);
        h = singfun(rsp, [exps(1) fracM-1]);
        y(k) = I + ((x(k)+1)/2)^fracM*sum(h);
    end
end

y = reshape(y, oldSize);

end