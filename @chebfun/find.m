function varargout = find(f)
%FIND   Find locations of nonzeros in a CHEBFUN.
%   FIND(F) returns a vector of all values at which the CHEBFUN F is nonzero.
%
%   [R, C] = FIND(F) returns two column vectors of the same length such that
%   [F(R(n),C(n)) for all n = 1:length(R)] is the list of all nonzero values of
%   the CHEBFUN F. One of the outputs holds dependent variable values, and
%   the other holds CHEBFUN row or column indices.
%
%   If the set of nonzero locations is not finite, an error is thrown.
%
% Example:
%   f = chebfun(@sin, [0. 2*pi]);
%   format long, find(f == 1/2) / pi
%
% See also ROOTS, EQ.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Is this method really required?

% Empty case:
if ( isempty(f) )
    varargout = {[], []};
    return
end

if ( (size(f.funs{1}, 2) > 1) && (nargout < 2) )
    error('CHEBFUN:CHEBFUN:find:arrout', ...
        'Use two output arguments for array-valued CHEBFUN objects.');
end

% Initialise:
x = [];
idx = [];

% Loop over columns:
for j = 1:size(f.funs{1}, 2)
    fj = extractColumns(f, j);

    % Loop over FUNs:
    for k = 1:numel(fj.funs)
        if ( ~iszero(fj.funs{k}) )
            % Continuous part is not identically zero!
            error('CHEBFUN:CHEBFUN:find:infset', ...
                'Nonzero locations are not a finite set.')
        end
    end
    
    xnew = fj.domain(fj.pointValues(:,1) ~= 0);   % Does all the real work.
    x = [x, xnew];                               %#ok<AGROW>
    idx = [idx, repmat(j, size(xnew))];          %#ok<AGROW>
end

if ( nargout == 1 )
    % Output has same shape as input.
    if ( ~f.isTransposed )
        x = x.';
    end
    varargout = {x};
else
    % Output is always column, but order matters.
    if ( ~f.isTransposed )
        varargout = {x.',idx.'};
    else
        varargout = {idx.',x.'};
    end
end

end
