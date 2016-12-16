function comet3(f, g, h, ignored)
%COMET3   Three-dimensional comet plot.
%   COMET3(F, G, H) displays a comet graph of the CHEBFUNs F, G, and H.
%
%   A comet graph is an animated graph in which a thick dot (the comet head)
%   traces the data points on the screen.
%
% See also COMET.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( (nargin == 1) && (numColumns(f) == 3) )
    % Case of a single quasimatrix input:
    f1 = extractColumns(f, 1);
    f2 = extractColumns(f, 2);
    f3 = extractColumns(f, 3);
    comet3(f1, f2, f3);
    return
end

if ( nargin > 3 )
    warning('CHEBFUN:CHEBFUN:comet3:nargin', ...
        'Fourth input to @CHEBFUN/COMET() is ignored.');
end

data = plotData(f, g, h);

comet3(data.xLine, data.yLine, data.zLine, 0);

end
