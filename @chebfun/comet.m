function comet(f, g, ignored)
%COMET   Two-dimensional comet plot.
%   COMET(F) displays a comet graph of the CHEBFUN F, and COMET(F, G) displays a
%   comet graph of the CHEBFUN F versus the CHEBFUN G.
%
%   A comet graph is an animated graph in which a thick dot (the comet head)
%   traces the data points on the screen.
%
% See also COMET3.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
%  For simplictiy, we simply extract the X and Y data from the standard PLOT()
%  method and then call the built-in COMET() method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a new figure
dummy = figure();
if ( nargin == 1 )
    % COMET(F)
    h = plot(f);
else
    if ( nargin > 2 )
        warning('CHEBFUN:comet:nargin', ...
            'Third input to @CHEBFUN/COMET() is ignored.');
    end
    % COMET(F, G)
    h = plot(f, g);
end

% Extract the data:
x = get(h, 'XData');
y = get(h, 'YData');

% Close the dummy figure:
close(dummy)

% Call the built-in COMET() method:
comet(x, y, 0);

end
