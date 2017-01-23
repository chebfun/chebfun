function out = ishappy(f)
%ISHAPPY   Happiness state of a CHEBFUN object.
%   ISHAPPY(F) returns TRUE of F is a happy CHEBFUN or quasimatrix and FALSE
%   otherwise.  F is happy if all of the FUNs that form each of its columns (or
%   rows, if F is a row CHEBFUN or quasimatrix) are happy.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

for j = 1:numel(f)
    for k = 1:numel(f(j).funs)
        out = all(get(f(j).funs{k}, 'ishappy'));
        if ( ~out )
            return
        end
    end
end

end
