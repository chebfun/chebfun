function ploteigenmodes(guifile, handles, selection, h1, h2)
% Plot the eigenmodes in the GUI

% TODO:  Documentation.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/chebfun/ for Chebfun information.

% selection == 0 corresponds to no selection being made, i.e. plot everything
if ( nargin < 3 )
    selection = 0;
end

if ( nargin < 4 )
    h1 = handles.fig_sol;
end

if ( nargin < 5 )
    h2 = handles.fig_norm;
end

% No recent solution available
if ( ~handles.hasSolution )
    return
end

% Obtain the most recent eigenvalues and eigenvectors
D = handles.latest.solution;
V = handles.latest.solutionT;

% Always create the same number of colours to preserve colours if selection
% is changed.
C = get(0, 'DefaultAxesColorOrder');
C = repmat(C, ceil(size(D)/size(C, 1)), 1);

% Number of unknown variables in the problem
numVar = size(V, 1);

% Number of columns in each CHEBMATRIX entry
numCol = size(V{1}, 2);


% Need to trim the data we are plotting if user has made a selection
if ( selection )
    % Pick out the selected eigenvalues
    D = D(selection);
    
    % Go through the rows of the CHEBMATRIX V and pick out the selected entries
    chebfunSelection = cell(numVar, 1);
    
    % Loop through the rows of the CHEBMATRIX V
    for selCounter = 1:numVar
        Vtemp = V{selCounter};
        chebfunSelection{selCounter} = Vtemp(:,selection);
    end
    
    % Convert the cell of selected CHEBFUNS to a CHEBMATRIX
    V = chebmatrix(chebfunSelection);

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

    if ( guifile.options.grid )
        grid on
    end

    title('Eigenvalues');
    xlabel('real');
    ylabel('imag');

    if ( any(selection) && (nargin < 4) )
        xlim(h1, xlim_sol);
        ylim(h1, ylim_sol);
    end

    axis equal
end

if ( isempty(h2) )
    return
end

isSystem = numVar > 1;

nV = numel(V);

realplot = get(handles.button_realplot, 'Value');
W = V;
if ( realplot )
    for k = 1:numVar
        V(k) = real(V{k});
    end
    s = 'Real part of eigenmodes';
else
    for k = 1:nV
        V(k) = imag(V{k});
    end
    s = 'Imaginary part of eigenmodes';
% TODO:  If this is really no longer necessary, remove this.
% else
%     V = V.*conj(V);
%     s = 'Absolute value of eigenmodes';
end

% TODO: This used to be a chebfunpref('plot_numpts'), do we still want to allow
% that?
maxPlotPoints = 2001;

axes(h2)
% set(h2,'NextPlot','add')
set(h2, 'ColorOrder', C)
if ( any(selection) && (nargin < 4) )
    xlim_norm = xlim(h2);
    ylim_norm = ylim(h2);
end

if ( ~isSystem )
    if ( (numCol == 1) && ~isreal(chebfun(W)) && ~isreal(1i*chebfun(W)) )
        xx = union(linspace(V.ends(1), V.ends(end), maxPlotPoints), V.ends);
        
        %TODO: Is there no slicker way to evaluate an array-valued CHEBFUN?
        Wxx = feval(chebfun(W), xx);
        
        WW = abs(reshape(Wxx, length(xx), size(V, 2)));
        
        plot(V(:,1), '-', 'linewidth', 2, 'color', C(1,:));
        hold on
        plot(xx, WW, '-', xx, -WW, '-', 'linewidth', 1, 'color', 'k');
        hold off
        xLims = V(:,1).domain;
    else
        % Convert to a CHEBFUN
        V = chebfun(V);
        for k = 1:size(V,2)
            plot(V(:,k), 'linewidth', 2, 'color', C(k,:));
            hold on
        end
        xLims = V(:,k).domain;
        hold off
    end

    if ( guifile.options.grid )
        grid on
    end

    ylabel(handles.varnames);
else
    LS = repmat({'-', '--', ':', '-.'}, 1, ceil(numVar/4));
    ylab = [];
    if ( (numCol == 1) && ~isreal(W{1}) && ~isreal(1i*W{1}) )
        V1 = V{1};
        xx = union(linspace(V1.ends(1), V1.ends(end), maxPlotPoints), V1.ends);
        for selCounter = 1:nV
            WW = abs(W{selCounter}(xx));
            plot(real(V{selCounter}), '-', 'linewidth', 2, 'linestyle', ...
                LS{selCounter});
            hold on
            plot(xx, WW, 'k', xx, -WW, 'k', 'linestyle', LS{selCounter});
        end
        xLims = V{selCounter}.domain;
        hold off
    else
        for selCounter = 1:numVar
            % If we are plotting selected e-funs, we need to pick out the colors
            if ( any(selection) )
                for sCounter = 1:length(selection)
                    plot(real(V{selCounter}(:,sCounter)), 'linewidth', 2, ...
                        'linestyle', LS{selCounter}, 'Color', C(sCounter,:));
                    hold on
                end
                xLims = V{selCounter}(:,sCounter).domain;
            else
                plot(real(V{selCounter}), 'linewidth', 2, 'linestyle',  ...
                    LS{selCounter});
                hold on
                xLims = V{selCounter}(:,1).domain;
            end
            ylab = [ylab handles.varnames{selCounter} ', ' ];
        end
    end
    hold off
    ylabel(ylab(1:end-2));
end
if ( any(selection) && (nargin < 4) )
    xlim(xlim_norm);
else
    Vdom = V.domain;
    xlim([Vdom(1) Vdom(end)]);
end
set(h2, 'NextPlot', 'replace')

xlabel(handles.indVarName);

% Set the xlim according to the domain of the function
title(s);
