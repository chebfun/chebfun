function varargout = plot(varargin)
%PLOT   Basic linear plot for DELTAFUN objects. 
%   PLOT(F) plots the DELTAFUN object F.
%
% See CHEBFUN/PLOT for details
%
% See also PLOTDATA, PLOT3.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% Deal with an empty input:
if ( isempty(varargin{1}) )
    if ( nargout == 1 )
        varargout{1} = plot([]);
    end
    return
end

% Store the hold state of the current axis:
holdState = ishold;

%%
% Plot the CHEBFUN:
f = varargin{1};
varargin(1) = [];
[h1, h2] = plot( f.funPart, varargin{:} );

if ( ~isempty(f.deltaLoc) && ~isempty(f.deltaMag) )
    deltaLoc = f.deltaLoc;
    deltaMag = f.deltaMag;
    hold on
    for i = 1:length(deltaLoc)
        plot( [deltaLoc(i), deltaLoc(i)], [0, deltaMag(1,i)], '-' );
        if ( deltaMag(1, i) > 0 )
            deltaMarker = '^';
        elseif ( deltaMag(1, i) < 0 )
            deltaMarker = 'v';
        end
        plot( deltaLoc(i), deltaMag(1, i), deltaMarker, 'MarkerFaceColor', 'b' );
    end
end
    
if ( holdState ) 
    hold on
else
    hold off
end

if ( nargout > 0 )
    varargout = {h1; h2};
end

end
