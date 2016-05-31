function h = openFigInCurrentFigure(figfile)
%OPENFIGINCURRENTFIGURE    Wrapper for OPENFIG to open in current figure.
%   H = OPENFIGINCURRENTFIGURE(FIGFILE) calls OPENFIG to open the *.fig
%   figure specified by the string FIGFILE but does some additional work
%   to force it to be opened in the current figure instead of in a new one.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Developer note:  This function exists to allow a user to force the
% special GUI-based plots for chebfun3 (e.g., isosurface(), slice(), etc.)
% to appear in a figure of choice (e.g., by calling figure() to set the
% current figure) instead of appearing in a new figure every time, just
% like ordinary MATLAB plots.

% Get a handle to the current figure and make it invisible while we set
% things up.
hgcf = gcf();
visibilityState = get(hgcf, 'Visible');
set(hgcf, 'Visible', 'off');

% Get rid of everything in the current figure but the menu and the
% toolbar if they exist.
for ( c = allchild(hgcf).' )
    switch ( class(c) )
        case 'matlab.ui.container.Menu'
        case 'matlab.ui.container.Toolbar'
        otherwise
            delete(c);
    end
end

% Open the *.fig file, again in an invisible window.
h = openfig(figfile, 'invisible');

% Copy the graphics elements from the *.fig figure into the current figure.
copyobj(allchild(h), hgcf);

% We don't need the *.fig figure anymore, so delete it.
delete(h);

% Return a handle to the current figure.
h = hgcf;

% Restore the visibility state.
set(h, 'Visible', visibilityState);

end