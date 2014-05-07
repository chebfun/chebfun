function plot(varargin)
%PLOT   Plot for CHEBMATRIX objects.
%   PLOT(A) plots the CHEBMATRIX object A.
%
%   If A contains only CHEBFUN and DOUBLE objects, A is converted to a
%   QUASIMATRIX, and CHEBFUN/PLOT is called. In this case:
%   
%   PLOT(A, S) allows various line types, plot symbols, and 
%   colors to be used, where S is a character string.
%
%   If A contains inf x inf blocks, CHEBMATRIX/SPY is called. In this case:
%    
%   SPY(A, DIM) uses the dimension vector DIM to create the picture.
%
%   SPY(A, DIM, DISCTYPE) uses a the chebDiscretization constructor DISCTYPE for
%   the visualization.
%
%   See also CHEBMATRIX, CHEBFUN/PLOT, CHEMATRIX/SPY.

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Deal with an empty input:
if ( isempty(varargin{1}) )
    if ( nargout == 1 )
        varargout{1} = plot([]);
    end
    return
end

A = varargin{1}; 
s = cellfun(@(b) min(size(b)), A.blocks);

% If A contains inf x inf blocks, call SPY().
if ( ~all(isfinite(s(:))) )
    spy(A, varargin{2:end});
% If A contains only CHEBFUN or DOUBLE, convert it to a QUASIMATRIX, and
% cal CHEBFUN/PLOT().
else
    A = chebfun(A);
    plot(A, varargin{2:end});
end

end