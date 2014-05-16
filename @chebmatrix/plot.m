function varargout = plot(A, varargin)
%PLOT   Plot for CHEBMATRIX objects.
%   PLOT(A) plots the CHEBMATRIX object A.
%
%   If A contains only CHEBFUN and DOUBLE objects, A is converted to a
%   QUASIMATRIX, and CHEBFUN/PLOT is called. In this case PLOT(A, S) allows
%   various line types, plot symbols, and colors to be used, where S is a
%   character string. See CHEBFUN/PLOT() for further details.
%
%   If A contains inf x inf blocks, CHEBMATRIX/SPY is called. In this case%
%   SPY(A, DIM, DISCTYPE) uses the dimension vector DIM and the discretization
%   DISCTYPE for the visualization. See CHEBMATRIX/SPY() for further details.
%
% See also CHEBFUN/PLOT, CHEMATRIX/SPY.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with an empty input:
if ( isempty(A) )
    if ( nargout == 1 )
        varargout{1} = plot([]);
    end
    return
end

s = cellfun(@(b) min(size(b)), A.blocks);
isQuasi = all(isfinite(s(:)));

if ( ~isQuasi )
    % If A contains inf x inf blocks, call SPY():
    h = spy(A, varargin{:});
    
else
    % If A contains only CHEBFUN or DOUBLE, convert it to a quasimatrix and
    % call CHEBFUN/PLOT():
    A = quasimatrix(A.blocks);
    h = plot(A, varargin{:});
    
end

if ( nargout > 0 )
    varargout = {h};
end

end