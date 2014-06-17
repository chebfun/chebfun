function varargout = chebpolyplot(A, varargin)
%CHEBPOLYPLOT   CHEBPOLYPLOT for CHEBMATRIX objects.
%   CHEBPOLYPLOT(A) plots the Chebyshev coefficients CHEBMATRIX object A. If A
%   contains only CHEBFUN and DOUBLE objects, A is converted to a QUASIMATRIX,
%   and CHEBFUN/CHEBPOLYPLOT() is called. In this case CHEBPOLYPLOT(A, S) allows
%   various line types, plot symbols, and colors to be used, where S is a
%   character string. See CHEBFUN/CHEBPOLYPLOT() for further details.
%
%   If A contains inf x inf blocks, an error is thrown.
%
% See also CHEBFUN/CHEBPOLYPLOT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with an empty input:
if ( isempty(A) )
    [varargout{1:nargout}] = plot([]);
    return
end

s = cellfun(@(b) min(size(b)), A.blocks);
isQuasi = all(isfinite(s(:)));

if ( ~isQuasi )
    % If A contains inf x inf blocks, call SPY():
    error('CHEBFUN:CHEBMATRIX:chebpolyplot:notQuasi', ...
        'CHEBPOLYPLOT does not support CHEBMATRIX objects of size INFxINF.');
end
     
% If A contains only CHEBFUN or DOUBLE, convert it to a quasimatrix and
% call CHEBFUN/CHEBPOLYPLOT():
A.blocks = reshape(A.blocks, 1, numel(A.blocks));
A = quasimatrix(A.blocks);
[varargout{1:nargout}] = chebpolyplot(A, varargin{:});

end
