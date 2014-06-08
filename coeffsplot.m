function varargout = coeffsplot(f, varargin)
%COEFFSPLOT   Display CHEBFUN coefficients graphically.
%   COEFFSPLOT(F) plots the underlying coefficients of the CHEBFUN F on a
%   semilogy scale. A horizontal line at the epslevel of F is also plotted. If
%   F is an array-valued CHEBFUN or has breakpoints, then a curve is plotted
%   for each FUN of each component (column) of F.
%
%   COEFFSPLOT(F, S) allows further plotting options, such as linestyle,
%   linecolor, etc, in the standard MATLAB manner. If S contains a string
%   'LOGLOG', the coefficients will be displayed on a log-log scale. If S
%   contains a string 'NOEPSLEVEL' the epslevel is not plotted.
%
%   H = COEFFSPLOT(F) returns a column vector of handles to lineseries
%   objects. The final entry is that of the epslevel plot.
%
% See also CHEBFUN/PLOT CHEBFUN/CHEBCOEFFSPLOT CHEBFUN/FOURCOEFFSPLOT

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with an empty input:
if ( isempty(f) )
    if ( nargout == 1 )
        varargout{1} = plot([]);
    end
    return
end

if isa(f.funs{1}.onefun,'chebtech')
    chebcoeffsplot(f,varargin{:});
elseif isa(f.funs{1}.onefun,'fourtech')
    fourcoeffsplot(f,varargin{:});
else
    error('CHEBFUN:coeffsplot:NotAvailable','Unknown chebfun type');
end

end

