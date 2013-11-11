function [x,w] = points(disc,kind)

d = disc.domain;
numint = disc.numIntervals;
n = disc.dimension;
if ( nargin < 2 )
    kind = 2;
end

if (numel(n) == 1)
    n = repmat(n,1,numint);
end

x = cell(numint,1);
w = cell(1,numint);
for k = 1:numint
    % Don't call again unless the dimension has changed.
    if ( k==1 ) || ( n(k) ~= n(k-1) )
        if ( kind == 2 )
            [x0,w0] = chebtech2.chebpts(n(k));
        else
            [x0,w0] = chebtech1.chebpts(n(k));
        end
    end
    dif = (d(k+1)-d(k))/2;
    x{k} = x0*dif + (d(k+1)+d(k))/2;
    if nargout > 1
        w{k} = w0*dif;
    end
end
x = cell2mat(x);
if (nargout > 1)
    w = cell2mat(w);
end
end
