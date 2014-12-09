function plotEigenmodes(handles, selection, h1, h2)
%PLOTEIGENMODE   Plot the eigenmodes in the GUI
% Calling sequence
%   PLOTEIGENMODES(HANDLES, SELECTION, H1, H2)
% where
%   HANDLES:    MATLAB handle object of the CHEBGUI figure.
%   SELECTION:  The user choice of a desired eigenvalue to be plotted, i.e. if 
%               a user selects an eigenvalue from the list shown after the
%               problem is solved.
%   H1:         A handle to the top plot of the CHEBGUI figure.
%   H2:         A handle to the bottom plot of the CHEBGUI figure.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% No recent solution available
if ( ~handles.hasSolution )
    return
else
    % Obtain the most recent eigenvalues and eigenvectors
    D = handles.latest.solution;
    V = handles.latest.solutionT;
end

% selection == 0 corresponds to no selection being made, i.e. plot everything
if ( nargin < 2 )
    selection = 0;
end

% Default figures to be plotted at.
if ( nargin < 3 )
    h1 = handles.fig_sol;
end
if ( nargin < 4 )
    h2 = handles.fig_norm;
end

% Always create the same number of colours to preserve colours if selection
% is changed.
C = get(0, 'DefaultAxesColorOrder');
C = repmat(C, ceil(length(D)/size(C, 1)), 1);

% Number of unknown variables in the problem
numVar = size(V, 1);

% Need to trim the data we are plotting if user has made a selection
if ( selection )
    % Pick out the selected eigenvalues
    D = D(selection);
    
    % Go through the rows of the CHEBMATRIX V and pick out the selected entries
    V = V(:,selection);

    % Pick out the colour needed for plotting the selected eigenvalues.
    C = C(selection,:);
end

if ( ~isempty(h1) )
    % Ensure that we still have the same x and y-limits on the plots. Only
    % do that when we are not plotting all the information

    if ( selection )
        xlim_sol = xlim(h1);
        ylim_sol = ylim(h1);
    end
    
    axes(h1)
    for k = 1:size(D)
        plot(real(D(k)), imag(D(k)), '.', 'markersize', 25, 'color', C(k,:));
        hold on
    end
    hold off

    % Show grid?
    if ( handles.guifile.options.grid )
        grid on
    end

    set(handles.panel_figSol, 'title', 'Eigenvalues (imag vs real parts)')
    
    set(h1, 'Fontsize', handles.fontsizePanels);

    if ( any(selection) && (nargin < 4) )
        xlim(h1, xlim_sol);
        ylim(h1, ylim_sol);
    end
    
end

if ( isempty(h2) )
    return
end

% Do we have a coupled system?
isSystem = numVar > 1;

% Do we want to plot the real or the imaginary parts of the eigenvalues?
realplot = get(handles.button_realplot, 'Value');
W = V;
if ( realplot )
    V = real(V);
    s = 'Real part of eigenmodes';
else
    V = imag(V);
    s = 'Imaginary part of eigenmodes';
end

% The number of points we use for plotting:
maxPlotPoints = 2001;

axes(h2)
set(h2, 'ColorOrder', C)
if ( any(selection) && (nargin < 4) )
    xlim_norm = xlim(h2);
end

% Do the plotting for the bottom figure. Coupled systems are more tricky than
% scalar problems.
if ( ~isSystem )
    % Deal with different kinds of plotting required depending on whether we
    % have real+imaginary parts or not.
    if ( (length(selection) == 1) && (selection > 0) && ~isreal(W{1}) && ~isreal(1i*W{1}) )
        d = V.domain;
        xx = union(linspace(d(1), d(end), maxPlotPoints), d).';
        WW = abs(feval(W{1}, xx));       
        plot(V{1}, '-', 'LineWidth', 2, 'color', C(1,:)); hold on
        plot(xx, WW, '-', xx, -WW, '-', 'LineWidth', 1, 'color', 'k'); hold off
    else
        for k = 1:size(V,2)
            plot(V{k}, 'LineWidth', 2, 'color', C(k,:));
            hold on
        end
        hold off
    end

    % Show grid?
    if ( handles.guifile.options.grid )
        grid on
    end
        
else
    % Linestyles for the eigenmodes.
    LS = repmat({'-', '--', ':', '-.'}, 1, ceil(numVar/4));
    % Label for the y-axis.
    ylab = [];
%     V = real(V);
    for varCounter = 1:numVar
        % If we are plotting selected e-funs, we need to pick out the colors
        if ( any(selection) )
            for sCounter = 1:length(selection)
                plot(V{varCounter,sCounter}, 'LineWidth', 2, ...
                    'LineStyle', LS{varCounter}, 'Color', C(sCounter,:));
                hold on
            end
        else
            plot(V(varCounter,:), 'LineWidth', 2, 'LineStyle', LS{varCounter});
            hold on
        end
        ylab = [ylab handles.varnames{varCounter} ', ' ]; %#ok<AGROW>
    end
    hold off
    
    % ylabel:
    ylabel(ylab(1:end-2));
    
end

% Set the xlim:
if ( any(selection) && (nargin < 4) )
    xlim(xlim_norm);
else
    Vdom = V.domain;
    xlim([Vdom(1) Vdom(end)]);
end
set(h2, 'NextPlot', 'replace')

% Set title of the panel
set(handles.panel_figNorm, 'title', s)
set(h2, 'Fontsize', handles.fontsizePanels);
end
