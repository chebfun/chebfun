function G = cell2quasi(F)
%CELL2QUASI   Convert a cell array of CHEBFUN objects to a quasimatrix.
%   CELL2QUASI(F) converts the cell array F of CHEBFUN objects in to a
%   quasimatrix.
%
% See also QUASIMATRIX, CHEB2CELL.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~iscell(F) )
    error('CHEBFUN:cell2quasi:notacell', ...
        'Input must be a cell array of scalar CHEBFUN objects.');
end

% Ensure that we have a scalar Chebfun:
isCheb = cellfun(@(f) isa(f, 'chebfun'), F);
chebIdx = find(isCheb);
if ( isempty(chebIdx) )
    error('CHEBFUN:cell2quasi:notacheb', ...
        'Input must be a cell array of scalar CHEBFUN objects.');
end

minSizeF = cellfun(@(f) min(size(f)), F);
if ( max(minSizeF) > 1 )
    error('CHEBFUN:cell2quasi:notascalar', ...
        'Input must be a cell array of scalar CHEBFUN objects.');
end

% Initialise G and make sure it's a CHEBFUN:
G = F{chebIdx(1)};

% Check the domains:
dom = G.domain;
for k = chebIdx
    pass = domainCheck(dom, F{k});
    if ( ~pass )
        error('CHEBFUN:cell2quasi:domaion', 'Inconsistent domains.');
    end
end

% Assign entries of F to columns of G:
for k = numel(F):-1:1
    G(k) = F{k};
end

end
