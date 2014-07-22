function plotTree(tree, varargin)
% PLOT Plot the AD tree of an anon (done recursively).

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.


% Starting values for plotting
startx = .5;
deltax = .25;

% Layout the tree. The returned tree will store information about where
% it's supposed to be plotted. xy and str are used to display data boxes
% when the plot is clicked.
tree = layoutNodes(tree, startx, deltax, tree.height+1, tree.height+1);
% 

% Create a figure, and plot the first node
fig = gcf;
plot(tree.x, tree.y, 'bo');
hold on

% Turn on datacursormode to make the plot clickable
h = datacursormode(fig);

% Set the update function, using the nested function below.
% set(h,'UpdateFcn',@textfun,'SnapToDataVertex','on');

% Update function for the clicking of the plot
function txt = textfun(obj,event_obj)
    % Display 'Time' and 'Amplitude'
    pos = get(event_obj,'Position');
    
    % Find index of text to be display
    findLoc = find(((pos(1) == xy(:,1)) & (pos(2) == xy(:,2))) == 1);  
    txt = str{findLoc};
    
    % Get rid of the hyperlink from the txt string. If there is a
    % hyperlink, we know it's preceded by a % sign
    perc_sign_loc = strfind(txt,'%');
    % and it finishes at the end of the first line
    if ~isempty(perc_sign_loc)
        first_newline = min(strfind(txt,sprintf('\n')));
        txt(perc_sign_loc:first_newline-1) = [];
    end
end

maxDiffOrder = tree.diffOrder;

plotTreePlot(tree, h, maxDiffOrder, varargin{:});

% datacursormode on

xlim([0 1])
ylim([0 1])
hold off
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Function evaluation tree')
end


function tree = layoutNodes(tree, treex, deltax, currHeight, maxHeight)


if ( ~isstruct(tree) )
    treeVal = tree;
    tree = [];
    tree.numArgs = 0;
    if ( isnumeric(treeVal) )
        tree.method = sprintf('%2.2f', treeVal);
    else
        tree.method = 'chebfun';
    end
    tree.y = (currHeight/(maxHeight+1));
    tree.x = treex;
    tree.diffOrder = 0;
    return
end

tree.y = (currHeight/(maxHeight+1));
tree.x = treex;
    
% Lay out nodes recursively
switch tree.numArgs
    case 0
        % Do nothing        
    case 1
        tree.center = layoutNodes(tree.center, treex, deltax, ...
            currHeight - 1, maxHeight);
    case 2
        tree.left  = layoutNodes(tree.left, treex - deltax, ...
            deltax/2, currHeight - 1, maxHeight);
        tree.right  = layoutNodes(tree.right, treex + deltax, ...
            deltax/2, currHeight - 1, maxHeight);
    case 3
        deltax = deltax / 2;
        [tree.left   xyl strl]  = layoutNodes(tree.left,treex-deltax,deltax/2,currHeight-1,maxHeight);
        [tree.center xyc strc] = layoutNodes(tree.center,treex,deltax/2,currHeight-1,maxHeight);
        [tree.right  xyr strr]  = layoutNodes(tree.right,treex+deltax,deltax/2,currHeight-1,maxHeight);
        xy = [xy;xyl;xyc;xyr];
        str = [str;strl;strc;strr];
end
end

function plotTreePlot(tree, h, maxDiffOrder, varargin)
% 
% for k = 1:3
%     if tree.found(k)
%         col{k} = 'r';
%     else
%         col{k} = 'b';
%     end
% end
colours = [160, 255, 0; 27, 244, 125; 0, 185, 140; 49, 105, 255; ...
    114, 109, 255; 137, 0, 211; 220, 44, 148; 255, 51, 91; ...
    252, 103, 1; 255, 145, 0; 255, 216, 0]/255;
numCols = size(colours, 1);
col = {'b', 'b', 'b'};
switch tree.numArgs
    case 0
        % Plot the constructor leaves in different colours.
        if ( strcmp(tree.method, 'constr') )
            varID = find(tree.ID == 1);
            colID = colours(mod(varID, numCols) + 1, :);
            plot(tree.x,tree.y,'o', 'color', colID);
            plot(tree.x,tree.y,'.', 'markersize', 20, 'color', colID);
        end
        % Do nothing
    case 1
        plot(tree.center.x,tree.center.y,'bo');
        plot([tree.x tree.center.x],[tree.y tree.center.y],'-','color',col{2},varargin{:})
        plotTreePlot(tree.center,h,varargin{:})
    case 2
        % Set up colors. The part of the tree that corresponds to the derivative
        % of the highest order gets plotted in red.
        if ( tree.left.diffOrder == maxDiffOrder )
            colorLeft = 'b';
        else
            colorLeft = 'b';
        end
        
        if ( tree.right.diffOrder == maxDiffOrder )
            colorRight = 'b';
        else
            colorRight = 'b';
        end
        
        % Plot left sub-tree.
        plot(tree.left.x,tree.left.y,'o', 'color', colorLeft);
        plot([tree.x tree.left.x], [tree.y tree.left.y], '-', ...
            'color', colorLeft, varargin{:})
        plotTreePlot(tree.left, h, maxDiffOrder, varargin{:})
        
        % Plot right sub-tree.
        plot(tree.right.x, tree.right.y, 'o', 'color', colorRight);
        plot([tree.x tree.right.x], [tree.y tree.right.y], '-', ...
            'color',colorRight, varargin{:})
        plotTreePlot(tree.right, h, maxDiffOrder, varargin{:})
        
        % Here we could do some adjustments if we're not using the full
        % width. Store minx, maxx?
    case 3
        plot(h,tree.left.x,tree.left.y,'bo');
        plot([tree.x tree.left.x],[tree.y tree.left.y],'-','color',col{1},varargin{:})
        plotTreePlot(tree.left,h,varargin{:})
        plot(h,tree.center.x,tree.center.y,'bo');
        plot([tree.x tree.center.x],[tree.y tree.center.y],'-','color',col{2},varargin{:})
        plotTreePlot(tree.center,h,varargin{:})
        plot(h,tree.right.x,tree.right.y,'bo');
        plot([tree.x tree.right.x],[tree.y tree.right.y],'-','color',col{3},varargin{:})
        plotTreePlot(tree.right,h,varargin{:})
end

% Add the method name to the plot:
if ( strcmp(tree.method, 'constr') )
    % If we're at a constructor leave, change the text slightly:
    if ( length(tree.ID) == 1)
        varString = 'u';    
    else
        varString = sprintf('u%i', find(tree.ID == 1));
    end
    text(tree.x + 0.02, tree.y - 0.01, varString, 'Interpreter', 'none')
else
    text(tree.x + 0.02, tree.y - 0.01, tree.method, 'Interpreter', 'none')
end

end


