function a = any(f, dim)
%ANY    True if any element of a chebfun is a nonzero number. ANY ignores
%entries that are NaN (Not a Number).
%
%   ANY(X, DIM), where X is a quasimatrix, works down the dimension DIM. If DIM
%   is the chebfun (continuous) dimension, then ANY returns a logical column
%   vector (or row) in which the Jth element is TRUE if any element of the Jth
%   column (or row) is nonzero. Otherwise, ANY returns a chebfun which takes the
%   value 1 wherever any of the columns (or rows) of X are nonzero, and zero
%   everywhere else.
%
% See also CHEBFUN/ALL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information

% ANY along the continuous dimension:
if ( isempty(f) )
    a = false;    % Empty ==> any == false.
    return
end

% Rotate dim if f is transposed:
if ( nargin == 2 && f.isTransposed )
    dim = mod(dim, 2) + 1; % Maps 1 to 2 and 2 to 1.
end

% [TODO]: Implement ANY() across rows of array-valued Chebfuns.
if ( nargin == 2 && dim ~= 1 )
    error('CHEBFUN:any:nargin', ...
        'ANY() along discrete direction not yet implemented.');
end

% Check the impulses:
if ( any(f.impulses(1,:)) )
    a = true;     % Non-trivial impulses ==> any == true.
    return
end

% Query each of the funs:
for k = 1:numel(f.funs)
    if ( any(f.funs{k}) )
        % If we find a non-empty fun, the chebfun is nonempty.
        a = true;
        return
    end
end

% If we've got to here, there are no non-trivial values:
a = false;

end