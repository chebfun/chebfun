function varargout = pivotplot(f,varargin)
%PIVOTPLOT(F) semilogy plot of the pivot values
%
% PIVOTPLOT(F) semilogy plot of the Gaussian elimination pivots taken
% during the construction of the chebfun2 F.
%
% H = PIVOTPLOT(F) returns a handle H to the figure.
%
% PIVOTPLOT(F,S) allows further plotting options, such as linestyle,
% linecolor, etc. If S contains a string 'LOGLOG', the psudeo sig will be
% displayed on a log-log scale.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

scaleTypeLog = false;
doWeHoldOn = ishold;

if ( isempty(f) ) 
   error('CHEBFUN2:PIVOTPLOT', 'Empty chebfun2 has no pivots to plot');
end 


% Parse input arguments
if ( nargin > 1 )
    for j = 1:length(varargin)
        if ( strcmpi(varargin{j}, 'loglog') )
            scaleTypeLog = true;
            varargin(j) = [];
            break
        end
    end
end

%% 
% Plot a semilogy or loglog 
plotopts = varargin;
UK = [ {1:length(pivots(f)), abs(pivots(f))}, plotopts ]; % store

if ( ~scaleTypeLog )
    h = semilogy( UK{:} );            % semilogy plot
else
    h = loglog( UK{:} );              % loglog plot
end

%%

if ( ~doWeHoldOn )
    hold off;  % hold off if we can. 
end

% output handle if appropriate
if ( nargout ~=0 )
    varargout = {h};
end

end