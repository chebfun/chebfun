function disp(f)
%DISP   Display a CHEBFUN object.
%
% See also DISPLAY.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% If the 'format loose' setting is enabled, we print additional linebreaks:
loose = strcmp(get(0, 'FormatSpacing'), 'loose');

s = '';

% Trivial case: (empty chebfun)
if ( isempty(f) )
    fprintf('   empty chebfun\n');
    if ( loose )
        fprintf('\n');
    end
    return
end

% Convert array-valued CHEBFUN to a quasimatrix for unified output.
if ( (numColumns(f) > 1) && (numel(f) == 1) )
    f = cheb2quasi(f);
end

% Transpose information:
if ( f(1).isTransposed )
    columnString = 'row';
else
    columnString = 'column';
end

if ( numel(f) == 1 )
    % Call the column display subfunction:
    fprintf(colDisp(f, columnString));
else
    for k = 1:numel(f)
        % Call the column display subfunction:
        fprintf(colDisp(f(k), [columnString, int2str(k)]));
        % Final line break:
        if ( loose )
            fprintf('\n');
        end
    end
end

if ( loose )
    fprintf('\n');
end

end

function s = colDisp(f, columnString)
% Initialise s for this column:
s = '   chebfun ';
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
s = [s, sprintf('\n       interval       length     endpoint values %s\n', extraItem)];
len = zeros(numFuns, 1);
for j = 1:numFuns
    len(j) = length(f.funs{j});

    if ( ~isreal(f.funs{j}) )
        % For complex-valued funs, we don't display the values.

        % Print information to screen:
        s = [s, sprintf('[%8.2g,%8.2g]   %6i     complex values %s\n', ...
            f.domain(j), f.domain(j+1), len(j), extraData{j})];

    else

        % Grab values at endpoints:
        endvals = [get(f.funs{j}, 'lval'), get(f.funs{j}, 'rval')];

        % Tweak the endpoint values: (This prevents -0 and +0)
        if ( ~any(isnan(endvals)) )
            endvals(~logical(abs(endvals))) = 0;
        end

        % Print information to screen:
        s = [s, sprintf('[%8.2g,%8.2g]   %6i  %8.2g %8.2g %s\n', ...
            f.domain(j), f.domain(j+1), len(j), endvals, extraData{j})];

    end
end

% Display vscale:
s = [s, sprintf('vertical scale = %3.2g ', vscale(f, 'sup'))];

% Display total length for piecewise chebfuns:
if ( numFuns > 1 )
    s = [s, sprintf('   Total length = %i', sum(len))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display for delta functions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out = get(f, 'deltas');
if ( ~isempty(out) )
    deltaLoc = out(1,:);
    deltaMag = out(2:end,:);
    s = [s, sprintf('\nDelta functions:\n')];
    m = size(deltaMag, 1);
    n = size(deltaMag, 2);
    for i = 1: m
        for j = 1:n
            s = [s, sprintf('%8.2g', deltaMag(i,j))];
        end
        s = [s, sprintf('\n')];
    end
    s = [s, sprintf( 'Locations:\n')];
    s = [s, sprintf('%8.2g', deltaLoc)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = [s, '\n'];

end
