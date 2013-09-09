function [x, w, v] = chebpts(n, dom, type)

% TODO: Document this file.
% TODO: Return W and X.

if ( nargin == 2 )
    if ( length(dom) == 1 )
        type = dom;
        dom = chebfun.pref('domain');
    else
        type = chebtech.pref('tech');
        type = str2double(type(end));
    end
end
if ( nargin == 1 )
    dom = chebfun.pref('domain');
    type = chebtech.pref('tech');
    type = str2double(type(end));
end

if ( length(n) == 1 && length(dom) > 1 )
    n = repmat(n, 1, length(dom) - 1);
elseif ( length(n) ~= length(dom) - 1 )
    error
end

f = feval(['chebtech', num2str(type)]);

if ( length(n) > 1 )
    x = cell(length(n), 1);
    for k = 1:numel(n)
        xk = f.chebpts(n(k));
        x{k} = scale(xk, dom(k:k+1));
    end
    x = cell2mat(x);
else
    x = f.chebpts(n);
    x = scale(x, dom);
end

end

function y = scale(x, dom)
% Scale the nodes.
if ( dom(1) == -1 && dom(2) == 1 )
    y = x;
    return
end
a = dom(1);
b = dom(2);
y = b*(x + 1)/2 + a*(1 - x)/2;
end