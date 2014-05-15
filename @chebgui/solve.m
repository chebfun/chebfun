function varargout = solve(guifile)
% SOLVE Called when a user hits calls the solve method for a chebgui object
% outside the GUI (i.e. SOLVE(CG), where CG is a CHEBGUI object).

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/chebfun/ for Chebfun information.

% Call different solver methods, depending on the type of the problem.
if ( strcmpi(guifile.type, 'bvp') )
    [varargout{1}, varargout{2}] = solveguibvp(guifile);
elseif ( strcmpi(guifile.type, 'pde') )
    [varargout{1}, varargout{2}] = solveguipde(guifile);
    
    % If only one output, need to switch around output arguments:
    if ( nargout < 2 )
        varargout{1} = varargout{2};
    end
else
    [varargout{1}, varargout{2}] = solveguieig(guifile);

    % If only one output, need to switch around output arguments:
    if ( nargout < 2 )
        varargout{1} = diag(varargout{2});
    end
end