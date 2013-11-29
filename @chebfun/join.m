function f = join(varargin)
%JOIN   Join together two CHEBFUN objects.
%   F = JOIN(F1, F2, ...) joins together the CHEBFUN objects F1, F2, ..., to
%   create a piecewise CHEBFUN F on a larger domain. F.domain(1) = F1.domain(1),
%   F.domain(2) = F.domain(1) + length(F1), F.domain(2) = F.domain(1) +
%   length(F1) + length(F2), and so on.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: This method needs a test.

% Check that all inputs have the same transposition state.
transStates = cell2mat(cellfun(@(f) logical(f.isTransposed), varargin, ...
    'UniformOutput', false).');
if ( ~(all(transStates) || all(~transStates)) )
    error('CHEBFUN:join:trans', ...
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
