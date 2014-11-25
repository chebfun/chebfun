function plotTree(tree)
%PLOTTREE   Plot a syntax tree (as stored in a TREEVAR object).
%   Calling sequence:
%      TREEVAR.PLOTTREE(TREE)
%   where the input is:
%      TREE:   A MATLAB struct, describing the syntax tree of a mathematical
%              expression.
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

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.


% Starting values for plotting
startx = .5;
deltax = .25;

% Layout the tree. The returned tree will store information about where
% it's supposed to be plotted.
tree = layoutNodes(tree, startx, deltax, tree.height + 1, tree.height + 1); 

% Create a figure, and plot the first node
plot(tree.x, tree.y, 'bo');
hold on

% Maximum differential orders that appear in the tree
maxDiffOrder = tree.diffOrder;

% Start the recursive calls to the actual plotting method:
plotTreePlot(tree, maxDiffOrder);

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

function plotTreePlot(tree, maxDiffOrder)
%PLOTREEPLOT   The actual plotting of the syntax tree.

% Specify 20 different colours, so that each variable appearing in the tree can
% be plotted with its own colour:
colours = [160, 255, 0; 27, 244, 125; 0, 185, 140; 49, 105, 255; ...
    114, 109, 255; 137, 0, 211; 220, 44, 148; 255, 51, 91; ...
    252, 103, 1; 255, 145, 0; 255, 216, 0]/255;

% Store the number of colours we use (so that we can take mod() later on):
numCols = size(colours, 1);

% Plot, depending on the number of arguments that the method has.
switch tree.numArgs
    case 0
        % Plot the constructor leaves in different colours.
        if ( strcmp(tree.method, 'constr') )
            varID = find(tree.ID == 1);
            colID = colours(mod(varID, numCols) + 1, :);
            plot(tree.x, tree.y, 'o', 'color', colID);
            plot(tree.x, tree.y, '.', 'markersize', 20, 'color', colID);
        end
        
    case 1
        % Plot syntax tree for univariate methods:
        plot(tree.center.x, tree.center.y, 'bo');
        plot([tree.x tree.center.x], [tree.y tree.center.y], '-')
        plotTreePlot(tree.center, maxDiffOrder)
    
    case 2
        % Plot syntax tree for bivariate methods:        
        % Plot left sub-tree.
        plot(tree.left.x, tree.left.y, 'o');
        plot([tree.x tree.left.x], [tree.y tree.left.y], '-')
        plotTreePlot(tree.left, maxDiffOrder)
        
        % Plot right sub-tree.
        plot(tree.right.x, tree.right.y, 'o');
        plot([tree.x tree.right.x], [tree.y tree.right.y], '-')
        plotTreePlot(tree.right, maxDiffOrder)
end

% Add the method name to the plot:
if ( strcmp(tree.method, 'constr') )
    % If we're at a constructor leaf, change the text slightly:
    if ( length(tree.ID) == 1)
        % Scalar case.
        varString = 'u';    
    else
        % Systems case, print u and the variable number:
        varString = sprintf('u%i', find(tree.ID == 1));
    end
    text(tree.x + 0.02, tree.y - 0.01, varString, 'Interpreter', 'none')
else
    text(tree.x + 0.02, tree.y - 0.01, tree.method, 'Interpreter', 'none')
end

end


