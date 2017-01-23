function disp(L, name)
%DISP   Print summary of CHEBMATRIX contents.
%
% See also DISPLAY, CHEBMATRIX/SPY.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[m, n] = size(L);
c = blockClasses(L);

loose = strcmp(get(0, 'FormatSpacing'), 'loose');

fprintf('   %i x %i chebmatrix of block types:\n', m, n)

if ( loose )
    fprintf('\n');
end

% Lite mode, either no Java, or Matlab is not running in desktop mode.
if ( ~usejava('jvm') || ~usejava('desktop') )
    disp(c);
    return
end

% We need to store the ANS variable in the workspace so that it can be made
% clickable:
if ( nargin < 2 )
    name = 'ans';
end
if ( strcmp(name, 'ans') )
    name = ['chebmatrix_' num2str(ceil(1e6*rand(1)), 6)];
    assignin('base', name, L);
end

% Print information about each block:
ml = cellfun(@length, c);
for j = 1:m
    fprintf('    ');
    for k = 1:n
        ws = repmat(' ', 1, ml(k) - length(c{j,k}) + 4);
        fprintf('''<a href="matlab: disp(%s{%d,%d})">%s</a>''%s', ...
            name, j, k, c{j,k}, ws);
    end
    fprintf('\n');
end

if ( loose )
    fprintf('\n');
end

end
