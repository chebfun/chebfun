function varargout = plotcoeffs(A, varargin)
%PLOTCOEFFS   PLOTCOEFFS for CHEBMATRIX objects.
%   PLOTCOEFFS(A) plots the Chebyshev coefficients CHEBMATRIX object A. If A
%   contains only CHEBFUN and DOUBLE objects, A is converted to a QUASIMATRIX,
%   and CHEBFUN/PLOTCOEFFS() is called. In this case PLOTCOEFFS(A, S) allows
%   various line types, plot symbols, and colors to be used, where S is a
%   character string. See CHEBFUN/PLOTCOEFFS() for further details.
%
%   If A contains inf x inf blocks, an error is thrown.
%
% See also CHEBFUN/PLOTCOEFFS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with an empty input:
if ( isempty(A) )
    [varargout{1:nargout}] = plot([]);
    return
end

s = cellfun(@(b) min(size(b)), A.blocks);
isQuasi = all(isfinite(s(:)));

if ( ~isQuasi )
    % If A contains inf x inf blocks throw an error:
    error('CHEBFUN:CHEBMATRIX:plotcoeffs:notQuasi', ...
        'PLOTCOEFFS does not support CHEBMATRIX objects of size INFxINF.');
end
     
% If A contains only CHEBFUN or DOUBLE, convert it to a quasimatrix and call
% CHEBFUN/PLOTCOEFFS():
A.blocks = reshape(A.blocks, 1, numel(A.blocks));
A = quasimatrix(A.blocks);
[varargout{1:nargout}] = plotcoeffs(A, varargin{:});

end
