function display(f)
%DISPLAY   Display a chebfun object.
%   DISPLAY(F) outputs important information about the chebfun F to the command
%   window, including its domain of definition, its length (number of sample
%   values used to represent it), and a summary of its values at its endpoints.
%   DISPLAY(F) is called automatically when the semicolon is not used at the end
%   of a statement that results in a CHEBFUN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% If the 'format loose' setting is enabled, we print additional linebreaks:
loose = strcmp(get(0, 'FormatSpacing'), 'loose');

% Print the variable name if possible:
if ( loose )
    fprintf('\n%s = \n\n', inputname(1))
else
    fprintf('%s = \n', inputname(1))
end

% Trivial case: (empty chebfun)
if ( isempty(f) )
    fprintf('   empty chebfun\n')
    return
end

% Transpose information:
if ( f.isTransposed )
    columnString = '   chebfun row';
else
    columnString = '   chebfun column';
end

% Number of pieces (i.e., funs) information:
numFuns = numel(f.funs);
if ( numFuns > 1 )
    disp([columnString, ' (' int2str(numFuns) ' smooth pieces)'])
else
    disp([columnString, ' (1 smooth piece)'])
end

% Extra information:
[extraItem, extraData] = dispInfo(f);

% Loop through each of the funs to display the following information:
fprintf('       interval       length   endpoint values %s \n', extraItem)
len = zeros(numFuns,1);

for j = 1:numFuns
    len(j) = length(f.funs{j});

    if ( min(size(f)) > 1 )
        % For array-valued funs, we don't display the values.

        % Print information to screen:
        fprintf('[%8.2g,%8.2g]   %6i    array-valued (%d pieces)\n', ...
            f.domain(j), f.domain(j+1), len(j), size(f.funs{j}, 2));

    elseif ( ~isreal(f.funs{j}) )
        % For complex-valued funs, we don't display the values.

        % Print information to screen:
        fprintf('[%8.2g,%8.2g]   %6i    complex values %s \n', ...
            f.domain(j), f.domain(j+1), len(j), extraData{j});

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
%         endvals(abs(endvals) < zeroTol & isinf(f.domain(j:j+1))) = 0;
        % Cheat zeros on bounded domains!:
        endvals(abs(endvals) < zeroTol) = 0;

        % Print information to screen:
        fprintf('[%8.2g,%8.2g]   %6i  %8.2g %8.2g %s \n', ...
            f.domain(j), f.domain(j+1), len(j), endvals, extraData{j});

    end
end

% Display epslevel:
fprintf('Epslevel = %3.3g, ', epslevel(f))
fprintf('Vertical scale = %.3g', vscale(f))

% Display total length for piecewise chebfuns:
if ( numFuns > 1 )
    fprintf(', Total length = %i.', sum(len))
else
    fprintf('.')
end

% Final line break:
if ( loose )
    fprintf('\n\n')
else
    fprintf('\n')
end

end
