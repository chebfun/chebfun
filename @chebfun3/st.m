function varargout = st(f)
%ST     SLICE-TUCKER decomposition of a CHEBFUN3.
%
%   [CORE, C, R, T] = ST(F) produces a core tensor CORE of size rank(F) and 
%   quasimatrices C, R, and T of size inf-by-length(F,1), inf-by-length(F,2) 
%   and inf-by-length(F,3), respectively such that 
%   f(x,y,z) = CORE x_1 C(x,:,:) x_2 R(:,y,:) x_3 T(:,:,z).
%
%   CORE = ST(F) returns the core tensor used in the construction of F.

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