function f = simplify(f, varargin)
%SIMPLIFY   Simplfy a chebfun.
%   [TODO]: Document this.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Simplfy each of the FUN objects:
for k = 1:numel(f.funs)
    f.funs{k} = simplify(f.funs{k}, varargin{:});
end

end
