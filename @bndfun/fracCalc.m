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
    f = funs{1};
    f.onefun = (diff(funs{1}.domain)/2)^fracM*(1/gamma(fracM))*...
        singfun(@(x) opOnePiece(funs{1}.onefun, x, fracM), [exps(1)+fracM exps(2)]);
else
    f = funs{end};
    f.onefun = (diff(funs{1}.domain)/2)^fracM*(1/gamma(fracM))*...
        singfun(@(x) opMultiPiece(funs, x, fracM), exps);
end

end

%%
function y = opMultiPiece(funs, x, fracM)

% Get the exponents of G:
g = funs{end}.onefun;
exps = g.exponents;
sp = g.smoothPart;

oldSize = size(x);
x = x(:);
l = length(x);
y = zeros(l, 1);

numFuns = numel(funs);

subDoms = zeros(numFuns, 2);
maps = cell(1, numFuns);

for j = 1:numFuns
    subDoms(j, :) = funs{j}.domain;
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

%%
function y = opOnePiece(f, x, M)

% Get the exponents of F:
exps = f.exponents;
sp = f.smoothPart;

oldSize = size(x);
x = x(:);
l = length(x);
y = zeros(l, 1);

for k = 1:l
    
    if ( x(k) == -1 )
        % When sampled at -1, the definite integral given above, i.e.
        % sum(integrand(f, x)) should be zero, which is equivalent to a zero
        % SINGFUN:
        y(k) = 0;
    elseif ( x(k) == 1 )
        % Call SINGFUN constructor:
        g = singfun(sp, [exps(1) exps(2)+M-1]);
        y(k) = sum(g);
    else
        % Call SINGFUN constructor:
        rsp = restrict(sp, [-1 x(k)]);
        g = singfun(rsp, [exps(1) M-1]);
        y(k) = ((x(k)+1)/2)^M*sum(g);
    end
end

y = reshape(y, oldSize);

end