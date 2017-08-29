function varargout = tucker(f)
%TUCKER   SLICE-TUCKER expansion of a CHEBFUN3 object.
%   [CORE, COLS, ROWS, TUBES] = TUCKER(F) returns the core tensor CORE and 
%   the three factor quasimatrices COLS, ROWS, and TUBES in the low-rank 
%   representation of a CHEBFUN3 object F. The factor quasimatrices are of 
%   size Inf-by-length(F, 1), Inf-by-length(F, 2) and Inf-by-length(F, 3), 
%   respectively and we have
%
%   F(x,y,z) = CORE x_1 COLS(x,:,:) x_2 ROWS(:,y,:) x_3 TUBES(:,:,z).
%
%   CORE = TUCKER(F) returns the core tensor used in the construction of F.
%
% See also CHEBFUN2/CDR.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    varargout = cell(1, nargout); 
    return
end

% Get the low rank representation for f. 
fCore = f.core;
fCols = f.cols; 
fRows = f.rows;
fTubes = f.tubes;

% Output:
if ( nargout <= 1 )
    varargout = {fCore};
else
    % ST decomposition
    varargout = {fCore, fCols, fRows, fTubes};
end

end