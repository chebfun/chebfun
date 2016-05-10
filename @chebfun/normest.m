function out = normest(f)
% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Check this and document.

for j = 1:numel(f)
    out = 0;
    for k = 1:numel(f(j).funs);
        out = out + normest(f(j).funs{k});
    end
end

end
