function y = subsref(f,s)
%SUBSREF Subscript reference for singfun.
%   F(X) evaluates the singfun F at the point(s) in X. This is an alias
%   for FEVAL(F,X).
%
%   See also SINGFUN/FEVAL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

switch s(1).type
    case '()'
        x = s(1).subs;
        y = feval(f,x{1});
    otherwise
        error('Chebfun:singfun:invalidSubsref',...
            'Invalid subscripting reference.')
end

end

