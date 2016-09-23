function plotTree(tree, varNames)
%PLOTTREE   Plot a syntax tree (as stored in a TREEVAR object).
%   Calling sequence:
%      TREEVAR.PLOTTREE(TREE, VARNAMES)
%   where the input is:
%      TREE:     A MATLAB struct, describing the syntax tree of a mathematical
%                expression.
%      VARNAMES: An optional cellstring, whose elements are the name of the
%                variables that appear in a problem.
%
%   Usually, this method is called from within the TREEVAR plot() method.
%
%   Example:
%      % First, define a TREEVAR and carry out some operations:
%      u = treeVar();
%      v = cos(u);
%      w = sin(diff(u));
%      t = v + w;
%      % The following are equivalent:
%      plot(t)
%      treeVar.plotTree(t.tree)
%
% See also TREEVAR.PLOT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Number of variables involved in the problem
numVars = length(tree.ID);

% Default variable names if none are passed:
if ( nargin < 2 )
    varNames = cell(1, numVars + 1);
    varNames{1} = 't';
    if ( numVars == 1 )
        % Simpler output for scalar case
        varNames{2} = 'u';
    else
        for varCounter = 1:numVars
            varNames{varCounter + 1} = sprintf('u%i', varCounter);
        end
    end
end

% Starting values for plotting
startx = .5;
deltax = .25;

% Layout the tree. The returned tree will store information about where
% it's supposed to be plotted.
tree = layoutNodes(tree, startx, deltax, tree.height + 1, tree.height + 1); 

% Create a figure, and plot the first node
plot(tree.x, tree.y, '.','markersize', 25);
hold on

% Maximum differential orders that appear in the tree
maxDiffOrder = tree.diffOrder;

% Start the recursive calls to the actual plotting method:
plotTreePlot(tree, maxDiffOrder, varNames);

% Some final massaging of the figure:
xlim([0, 1])
ylim([0, 1])
hold off
set(gca, 'xtick', [])
set(gca, 'ytick', [])
title('Function evaluation tree')
end


function tree = layoutNodes(tree, treex, deltax, currHeight, maxHeight)
%LAYOUTNODES   Return coordinates for where to plot nodes of syntax trees.

if ( ~isstruct(tree) )
    % In case of scalars or CHEBFUNS appearing in the expression, the
    % corresponding nodes will not be syntax trees. Need to treat this case
    % separately.
    
    % Store the value of the current "tree":
    treeVal = tree;
    
    % Create a dummy tree, needed for storing the information required when the
    % tree is plotted:
    tree = [];
    tree.numArgs = 0;
    if ( isnumeric(treeVal) )
        tree.method = num2str(treeVal, '%2.2f');
    else
        tree.method = 'chebfun';
    end
    
    % Store coordinates where the node should be plotted:
    tree.y = currHeight/(maxHeight+1);
    tree.x = treex;
    
    % Scalars and CHEBFUNS have zero diffOrders.
    tree.diffOrder = 0;
    return
end

% Store coordinates where the current node should be plotted:
tree.y = currHeight/(maxHeight+1);
tree.x = treex;
    
% Lay out nodes recursively
switch tree.numArgs
    case 0
        % Do nothing        
    case 1
        % Current tree has one child:
        tree.center = layoutNodes(tree.center, treex, deltax, ...
            currHeight - 1, maxHeight);
    case 2
        % Current tree has two childs:
        tree.left  = layoutNodes(tree.left, treex - deltax, ...
            deltax/2, currHeight - 1, maxHeight);
        tree.right  = layoutNodes(tree.right, treex + deltax, ...
            deltax/2, currHeight - 1, maxHeight);
end

end

function plotTreePlot(tree, maxDiffOrder, varNames)
%PLOTREEPLOT   The actual plotting of the syntax tree.

% Nice modern Matlab colours for lines and nodes:
lineColor = [0 0.447 0.741];
depVarColor = [0.85 0.325 0.098];
indVarColor = [0.929 0.694 0.125];

% Plot, depending on the number of arguments that the method has.
switch tree.numArgs
    case 0
        % Plot the constructor leaves in different colours.
        if ( strcmp(tree.method, 'constr') )
            if ( ~any(tree.ID) )
                % Independent variable
                col = indVarColor;
            else
                col = depVarColor;
            end
            plot(tree.x, tree.y, '.','markersize', 25, 'color', col);
        end
        
    case 1
        % Plot syntax tree for univariate methods:
        plot(tree.center.x, tree.center.y, '.','markersize', 25, ...
            'color', lineColor);
        plot([tree.x tree.center.x], [tree.y tree.center.y], 'color',lineColor)
        plotTreePlot(tree.center, maxDiffOrder, varNames)
    
    case 2
        % Plot syntax tree for bivariate methods. Start with left sub-tree.
        plot(tree.left.x, tree.left.y, '.','markersize', 25, ...
            'color', lineColor);
        plot([tree.x tree.left.x], [tree.y tree.left.y], 'color',lineColor)
        plotTreePlot(tree.left, maxDiffOrder, varNames)
        
        % Plot right sub-tree.
        plot(tree.right.x, tree.right.y, '.','markersize', 25, ...
            'color', lineColor);
        plot([tree.x tree.right.x], [tree.y tree.right.y], 'color',lineColor)
        plotTreePlot(tree.right, maxDiffOrder, varNames)
end

% Add the method name to the plot:
if ( strcmp(tree.method, 'constr') )
    % If we're at a constructor leaf, change the text slightly:
    varID = find(tree.ID == 1);
    if ( isempty(varID) )
        % Independent variable
        varString = varNames{1};
    else
        varString = varNames{varID + 1};
    end
    
    text(tree.x + 0.02, tree.y - 0.01, varString, 'Interpreter', 'none', ...
        'fontsize',14)
else
    text(tree.x + 0.02, tree.y - 0.01, tree.method, 'Interpreter', 'none', ...
        'fontsize',14)
end

end


