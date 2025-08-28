function f = join(varargin)
%JOIN   Join together two or more CHEBFUN objects.
%   F = JOIN(F1, F2, ...) joins together the CHEBFUN objects F1, F2, ..., to
%   create a piecewise CHEBFUN F on a larger domain. F1, F2, ... must all have
%   the same transposition state; the output F will have the same transposition
%   state as the inputs. The left endpoint of the domain of F is F1.domain(1),
%   and the remaining points in the domain are determined by the adding on the
%   lengths of the successive subintervals forming the domains of F1, F2, etc.
%   For example, if F1 has domain [-1 -0.5 0] and F2 has domain [1 1.25 2],
%   then the domain of F will be [-1 -0.5 0 0.25 1].
%
%   The number of columns (or rows) in F1, F2, ... must be the same, else an
%   error is thrown.
%
% See also HORZCAT, VERTCAT, NEWDOMAIN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Trivial cases:
if ( nargin == 0 )
    f = chebfun();
    return
elseif ( nargin == 1 )
    f = varargin{1};
    return
end

% Deal with quasimatrices:
numCols = cellfun(@numColumns, varargin);
numEls = cellfun(@numel, varargin);
if ( any(numCols > 1) )
    args = cellfun(@cheb2cell, varargin, 'UniformOut', false);
    try
        args = reshape([args{:}], max(size(args{1})), nargin);
    catch
        error('CHEBFUN:CHEBFUN:join:dim', 'Matrix dimensions must agree.');
    end
    for k = numCols(1):-1:1
        f(k) = columnJoin(args{k,:});
    end
    if ( all(numEls == 1) )
        f = quasi2cheb(f);
    end
else
    f = columnJoin(varargin{:});
end

end

function f = columnJoin(varargin)

% Check that all inputs have the same transposition state.
transStates = cell2mat(cellfun(@(f) logical(f.isTransposed), varargin, ...
    'UniformOutput', false).');
if ( ~(all(transStates) || all(~transStates)) )
    error('CHEBFUN:CHEBFUN:join:columnJoin:trans', ...
        'All inputs to JOIN must have the same transposition state.');
end

% Remember the transposition state.
transState = transStates(1);

% Extract each of the FUNs:
funs = cellfun(@(f) f.funs, varargin, 'UniformOutput', false);
funs = [funs{:}];

% Extract the old domains and their lengths:
oldDom = cell2mat(cellfun(@(f) f.domain, funs, 'UniformOutput', false).');
lengths = diff(oldDom, [], 2);
newDom = oldDom;

% Update the domains on each of the FUNS:
for k = 2:numel(funs)
    newDom(k,:) = newDom(k-1,2) + [0, lengths(k)];
    funs{k} = changeMap(funs{k}, newDom(k,:));
end

% Combine in a CHEBFUN:
f = chebfun(funs);
f.isTransposed = transState;

end
