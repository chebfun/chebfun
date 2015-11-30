function F = quasi2cheb(F)
%QUASI2CHEB   Convert a quasimatrix to an array-valued CHEBFUN.
%   QUASI2CHEB(F) converts the quasimatrix F to an array-valued CHEBFUN by
%   taking the union of the domains of each of the columns:
%
% See also QUASIMATRIX, CHEB2QUASI, NUM2CELL.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( numel(F) < 2 )
    % F is already be a quasimatrix!
    return
end

% Unify the breakpoints:
F = restrict(F, get(F, 'domain'));

% Collect each column in a cell:
F = cheb2cell(F);

% Call CAT to collate the columns:
F = cat(2 - F{1}.isTransposed, F{:});

end
