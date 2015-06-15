function comet(f, g, ignored)
%COMET   Two-dimensional comet plot.
%   COMET(F) displays a comet graph of the CHEBFUN F, and COMET(F, G) displays a
%   comet graph of the CHEBFUN F versus the CHEBFUN G.
%
%   A comet graph is an animated graph in which a thick dot (the comet head)
%   traces the data points on the screen.
%
% See also COMET3.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 ) % COMET(F)
    data = plotData(f);
else               % COMET(F, G)
    if ( nargin > 2 )
        warning('CHEBFUN:CHEBFUN:comet:nargin', ...
            'Third input to @CHEBFUN/COMET() is ignored.');
    end
    data = plotData(f, g);
end

comet(data.xLine, data.yLine, 0);

end
