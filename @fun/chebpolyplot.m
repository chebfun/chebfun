function varargout = chebpolyplot(f, varargin)
%CHEBPOLYPLOT    Display Chebyshev coefficients graphically.
%
%   CHEBPOLYPLOT(F) assumes that the ONEFUN of F is based on Chebyshev
%   technology. If not, the method is expected to throw an error. The method
%   plots the Chebyshev coefficients of the ONEFUN of a FUN F on a semilogy
%   scale. A horizontal line at the EPSLEVEL of F.ONEFUN is also plotted. If F
%   is an array-valued FUN then a curve is plotted for each component (column) 
%   of F.
%
%   CHEBPOLYPLOT(F, S) allows further plotting options, such as linestyle,
%   linecolor, etc, in the standard MATLAB manner. If S contains a string
%   'LOGLOG', the coefficients will be displayed on a log-log scale. If S
%   contains a string 'NOEPSLEVEL', the EPSLEVEL is not plotted.
%
%   H = CHEBPOLYPLOT(F) returns a column vector of handles to lineseries
%   objects. The final entry is that of the EPSLEVEL plot.
%
%   This method is mainly aimed at Chebfun developers.
%
% See also CHEBPOLY, PLOT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with an empty input:
if ( isempty(f) )
    if ( nargout == 1 )
        varargout{1} = plot([]);
    end
    return
end


% Call CHEBPOLYPLOT() on the ONEFUN of F.
h = chebpolyplot(f.onefun, varargin{:});

% Give an output if one was requested:
if ( nargout > 0 )
    varargout{1} = h;
end

end