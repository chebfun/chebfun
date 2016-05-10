function varargout = fill(varargin)
%FILL  Filled 2-D CHEBFUN plots.
%   FILL(F, G, C) fills the 2-D region defined by CHEBFUN objects F and G with
%   the color specified by C. If necessary, the region is closed by connecting
%   the first and last point of the curve defined by F and G.
%
%   If C is a single character string chosen from the list 'r', 'g', 'b', 'c',
%   'm', 'y', 'w', 'k', or an RGB row vector triple, [r g b], the polygon is
%   filled with the constant specified color.
%
%   If F and G are array-valued CHEBFUN objects of the same size, one region per
%   column is drawn. FILL(F1, G1, C1, F2, G2, C2, ...) is another way of
%   specifying multiple filled areas. Note that FILL does not support
%   quasimatrix input.
%
%   FILL sets the PATCH object FaceColor property to 'flat', 'interp', or a
%   colorspec depending upon the value of the C matrix.
%
%   H = FILL(...) returns a column vector of handles to PATCH objects, one
%   handle per patch. The F, G, C triples can be followed by parameter/value
%   pairs to specify additional properties of the patches.
%
% See also AREA, PLOT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Quasimatrix support?

k = 1;
% Get discrete data from PLOTDATA():
while ( k < (length(varargin) - 1) )
    if ( isa(varargin{k}, 'chebfun') )
        if ( numel(varargin{k}) > 1 || numel(varargin{k+1}) > 1 )
            error('CHEBFUN:CHEBFUN:fill:quasi', ...
                'FILL does not support quasimatrices.');
        end
        % Call plotData():
        data = plotData(varargin{k}, varargin{k+1});
        % Remove NaNs (arrising from interior breakpoints):
        idx = ~any(isnan(data.xLine), 2);
        varargin(k:k+1) = {data.xLine(idx,:), data.yLine(idx,:)};
        k = k + 1;
    end
    k = k + 1;
end

if ( any(cellfun(@(f) isa(f, 'chebfun'), varargin)) )
    error('CHEBFUN:CHEBFUN:fill:oops', 'Unrecognised input sequence.');
end

% Call the built in FILL():
h = fill(varargin{:});

% Output handle:
if ( nargout == 1 )
    varargout{1} = h;
end

end

  
