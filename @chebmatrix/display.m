function display(L)
%DISPLAY   Print summary of CHEBMATRIX contents.
%   DISPLAY(L) prints the size of the CHEBMATRIX L and a list of the block class
%   types. If java is enabled, the class types for each block are hyperlinked,
%   and clicking them will call the display method for the item in that block.
%
% See also CHEBMATRIX.SPY.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

name = inputname(1);
[m, n] = size(L);

c = blockClasses(L);

loose = strcmp(get(0, 'FormatSpacing'), 'loose');
if ( loose )
    fprintf('\n');
end

fprintf('%s = %i x %i chebmatrix of block types:\n', name, m, n)

if ( loose )
    fprintf('\n');
end

% Lite mode, either no Java, or Matlab is not running in desktop mode.
if ( ~usejava('jvm') || ~usejava('desktop') )
    disp(c);
    return
end

% We need to store the ANS variable in the workspace, so that it can be made
% clickable:
if ( strcmp(name, 'ans') )
    name = ['chebmatrix_' num2str(ceil(1e6*rand(1)), 6)];
end

% Print information about each block>
ml = cellfun(@length, c);
for j = 1:m
    fprintf('    ');
    for k = 1:n
        ws = repmat(' ', 1, ml(k) - length(c{j,k}) + 4);
        fprintf('''<a href="matlab: display(%s{%d,%d})">%s</a>''%s', ...
            name, j, k, c{j,k}, ws);
    end
    fprintf('\n');
end

if ( loose )
    fprintf('\n');
end

end
