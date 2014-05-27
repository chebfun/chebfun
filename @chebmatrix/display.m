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

if ( ~usejava('jvm') || ~usejava('desktop') )
    disp(c);
    return
end

if ( strcmp(name, 'ans') )
    name = ['chabmatrix_' num2str(1e6*rand(1), 6)];
    assignin('base', name, L);
end

for j = 1:m
    for k = 1:n
        fprintf('      ''<a href="matlab: display(%s{%d,%d})">%s</a>''', ...
            name, j, k, c{j,k});
    end
    fprintf('\n');
end

if ( loose )
    fprintf('\n');
end

end
