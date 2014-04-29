function f = cheb2cell(f)
%CHEB2CELL  Convert columns of a quasimatrix or array-valued CHEBFUN to a cell.
%   C = CHEB2CELL(F) converts an INFxM array-valued CHEBFUN or quasimatrix F
%   into a 1xM cell array C by placing each column of F into a separate cell in
%   C. If F is an MxINF row CHEBFUN, then C is Mx1.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note:
%   This method is simply a wrapper for num2cell with a more descriptive name.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call NUM2CELL:
f = num2cell(f);

end
