function s = disp(f, name)
%DISPLAY   Display a chebfun object.
%   DISPLAY(F) outputs important information about the chebfun F to the command
%   window, including its domain of definition, its length (number of sample
%   values used to represent it), and a summary of its values at its endpoints.
%   DISPLAY(F) is called automatically when the semicolon is not used at the end
%   of a statement that results in a CHEBFUN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% If the 'format loose' setting is enabled, we print additional linebreaks:
loose = strcmp(get(0, 'FormatSpacing'), 'loose');

if ( nargin < 2 ) 
    name = inputname(1);
end

% Print the variable name if possible:
if ( loose )
    s = sprintf('\n%s = \n\n', name);
else
    s = sprintf('%s = \n', name);
end

% Trivial case: (empty chebfun)
if ( isempty(f) )
    s = [s, sprintf('   empty chebfun')];
    return
end

if ( numel(f) > 1 )
    for k = 1:numel(f)
        % Transpose information:
        if ( f(1).isTransposed )
            columnString = ['   chebfun row', int2str(k)];
        else
            columnString = ['   chebfun column', int2str(k)];
        end
        % Call the column display subfunction:
        s = [s, colDisp(f(k), columnString)];
        if ( k ~= numel(f) )
            s = [s, sprintf('\n')];
        end
        % Final line break:
        if ( loose )
            s = [s, sprintf('\n')];
        end
    end
    
elseif ( numColumns(f) > 1 )
    % Convert to a quasimatrx for unified output.
    f = cheb2quasi(f);
    s = disp(f, name);
    % TODO: Remove this!
    s = [s, sprintf('\n(P.S. I am an array-valued CHEBFUN!)')];
    return
    
else
    % Transpose information:
    if ( f.isTransposed )
        columnString = '   chebfun row';
    else
        columnString = '   chebfun column';
    end
    % Call the column display subfunction:
    s = [s, colDisp(f, columnString)];
end

if ( loose )
    s = [s, sprintf('\n')];
end

end

function s = colDisp(f, columnString)
% Initialise s for this column:
s = '';
% Number of pieces (i.e., funs) information:
numFuns = numel(f.funs);
if ( numFuns > 1 )
    s = [s, columnString, ' (' int2str(numFuns) ' smooth pieces)'];
else
    s = [s, columnString, ' (1 smooth piece)'];
end

% Extra information:
[extraItem, extraData] = dispData(f);

% Loop through each of the funs to display the following information:
s = [s, sprintf('\n       interval       length   endpoint values %s\n', extraItem)];
len = zeros(numFuns, 1);
for j = 1:numFuns
    len(j) = length(f.funs{j});

    if ( numColumns(f) > 1 )
        % For array-valued funs, we don't display the values.

        % Print information to screen:
        s = [s, sprintf('[%8.2g,%8.2g]   %6i    array-valued (%d pieces)\n', ...
            f.domain(j), f.domain(j+1), len(j), size(f.funs{j}, 2))];

    elseif ( ~isreal(f.funs{j}) )
        % For complex-valued funs, we don't display the values.

        % Print information to screen:
        s = [s, sprintf('[%8.2g,%8.2g]   %6i    complex values %s\n', ...
            f.domain(j), f.domain(j+1), len(j), extraData{j})];

    else

        % Grab values at endpoints:
        endvals = [get(f.funs{j}, 'lval'), get(f.funs{j}, 'rval')];

        % Tweak the endpoint values: (This prevents -0 and +0)
        if ( ~any(isnan(endvals)) )
            endvals(~logical(abs(endvals))) = 0;
        end

        % Cheat zeros:
        zeroTol = 100*max(get(f.funs{j}, 'vscale')*get(f.funs{j}, 'epslevel'));
        % Cheat zeros on unbounded domains:
%         endvals(abs(endvals) < zeroTol & isinf(f.domain(j:j+1))) = 0; %TODO
        % Cheat zeros on bounded domains!:
        endvals(abs(endvals) < zeroTol) = 0;

        % Print information to screen:
        s = [s, sprintf('[%8.2g,%8.2g]   %6i  %8.2g %8.2g %s\n', ...
            f.domain(j), f.domain(j+1), len(j), endvals, extraData{j})];

    end
end

% Display epslevel:
s = [s, sprintf('Epslevel = %i.', epslevel(f))];
s = [s, sprintf('  Vscale = %i.', vscale(f))];

% Display total length for piecewise chebfuns:
if ( numFuns > 1 )
    s = [s, sprintf('  Total length = %i.', sum(len))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dispaly for delta functions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[deltaMag, deltaLoc] = getDeltaFunctions(f);
if ( ~isempty(deltaMag) )    
    s = [s, sprintf('\nDelta functions:\n')];
    m = size(deltaMag, 1);
    n = size(deltaMag, 2);
    for i = 1: m
        for j = 1:n            
            s = [s, sprintf('%8.2g', deltaMag(i, j))];
        end
        s = [s, sprintf('\n')];
    end
    s = [s, sprintf( 'Locations:\n')];
    s = [s, sprintf('%8.2g', deltaLoc)];
    s = [s, sprintf('\n')];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
