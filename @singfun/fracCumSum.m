function f = fracCumSum(oneFuns, maps ,fracM)
%FRACCUMSUM   Fractional indefinite integral of a ONEFUN.
%  FRACCSUMSUM(ONEFUNS, MAPS, FRACM) returns the fractional indefinite integral
%  of order FRACM defined by the ONEFUNs stored in the cell ONEFUNS and their 
%  corresponding map structure which maps between the real domain and [-1 1]. 
%  Note that FRACM needs to be non-integer (fractional).
%
% See also CUMSUM, DIFF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Grab the number of ONEFUNs stored in ONEFUNS:
numOneFuns = numel(oneFuns);

% Check if each of ONEFUNS is a SINGFUN, if not type-cast it to a SINGFUN:
for k = 1:numOneFuns
    if ( ~isa(oneFuns{k}, 'singfun') )
        oneFuns{k} = singfun(oneFuns{k});
    end
end

% Grab the exponent of the last piece:
exps = oneFuns{end}.exponents;

% Exponents for the result:
newExps = exps;

% If there is only one piece, modify the exponents:
if ( numOneFuns == 1 )
    newExps(1) = newExps(1) + fracM;
end

% Call the SINGFUN constructor to construct the Riemann–Liouville fractional 
% integral:

f = (1/gamma(fracM))*singfun(@(x) op(oneFuns, maps, x, fracM), newExps);

end

%% The operator of the Riemann–Liouville fractional integral:
function y = op(oneFuns, maps, x, fracM)

% Get the smooth part and the exponents of the last piece:
g = oneFuns{end};
exps = g.exponents;
sp = g.smoothPart;

% Preallocation and vectorization:
oldSize = size(x);
x = x(:);
l = length(x);
y = zeros(l, 1);

numOneFuns = numel(oneFuns);

% Loop over each sample point:
for k = 1:l
        
    % Compute the contribution from the previous pieces. If there is only one 
    % piece, the contribution is zero:
    I = 0;
    
    for j = 1:(numOneFuns-1)
        
        t = maps{j}.inv(maps{end}.for(x(k)));
        kernel = singfun(@(s) (t-s).^(fracM-1), [0 0]);
        Ij = sum(kernel.*oneFuns{j});
        I = I + Ij;
    end
    
    if ( x(k) == -1 )
        % When sampled at -1, namely the left endpoint of the entire CHEBFUN or 
        % the left edge of a piece, set the value of the fractional integral 
        % at -1 as NaN and this value will be extrapolated during construction:
        y(k) = NaN;
        
    elseif ( x(k) == 1 )
        % Call SINGFUN constructor:
        h = singfun(sp, [exps(1) exps(2)+fracM-1]);
        y(k) = I + sum(h);
    else
        % Call SINGFUN constructor:
        rsp = restrict(sp, [-1 x(k)]);
        h = singfun(rsp, [exps(1) fracM-1]);
        y(k) = I + ((x(k)+1)/2)^(fracM + exps(1))*sum(h);
    end
    
end

% Reshape the result:
y = reshape(y, oldSize);

end