function [F, isArrayValued] = quasi2cheb(F)
%QUASI2CHEB   Convert a quasimatrix to an array-valued CHEBFUN.
%   QUASI2CHEB(F) converts the quasimatrix F to an array-valued CHEBFUN by
%   taking the union of the domains of each of the columns.
%
%   [G, ISARRAYVALUED] = QUASI2CHEB(F) returns a second output with the value
%   TRUE if G is a CHEBFUN (i.e,. numel(G) == 1), and FALSE if G remains a
%   quasimatrix.
%
% See also QUASIMATRIX, CHEB2QUASI, NUM2CELL.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( numel(F) < 2 )
    % F is already a quasimatrix!
    isArrayValued = 1;
    return
end

% Unify the breakpoints:
F = restrict(F, get(F, 'domain'));

% Collect each column in a cell:
F = cheb2cell(F);

% Call CAT to collate the columns:
F = cat(2 - F{1}.isTransposed, F{:});

isArrayValued = numel(F) == 1;

end
