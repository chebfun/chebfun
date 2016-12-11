function disp(F)
%DISP   Display a CHEBFUN3 to the command line.
% 
% See also CHEBFUN3/DISPLAY.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

loose = strcmp( get(0, 'FormatSpacing'), 'loose');

% Get display style and remove trivial empty CHEBFUN3 case.
if ( isempty(F) )
    fprintf('    empty chebfun3\n')
    if ( loose )
        fprintf('\n');
    end
    return
end

% Get information that we want to display:
dom = F.domain;                           % Domain
r = size(F.core); % For the zero function we still expect to 
% see that the core tensor is 1 x 1 x 1, rather than 0 x 0 x 0.
if ( numel(r) < 3 )
    r = [r, 1];
end
r1 = r(1); 
r2 = r(2); 
r3 = r(3);
vscl = vscale(F);                         % vertical scale
[m, n, p] = length(F);

% Check underlying tech is a TRIGTECH: 
techCol = get(F.cols.funs{1}, 'tech');
techRow = get(F.rows.funs{1}, 'tech');
techTube = get(F.tubes.funs{1}, 'tech');

% Display the information: 
if ( isa(techCol(), 'trigtech') && isa(techRow(), 'trigtech') ...
        && isa(techTube(), 'trigtech') )
    disp('   chebfun3 object  (trig)')
else
    disp('   chebfun3 object ')
end

% To align the display nicely, check how many digits the maximum rank has.
wid = num2str(max([length(num2str(r1)), length(num2str(r2)), ...
    length(num2str(r3))]));

if ( all(floor(dom) == dom) ) 
    % Corners of the domain are all integers.
    domainFormatString = 'domain: [%-d, %d] x [%-d, %d] x [%-d, %d]\n';
    
else    
    domainFormatString = ['domain: [%-.3g, %.3g] x [%-.3g, %.3g] x '...
        '[%-.3g, %.3g]\n'];
end

str = ['   cols: Inf x %-' wid 'd chebfun\n   rows: Inf x %-' wid 'd '...
    'chebfun\n  tubes: Inf x %-' wid 'd chebfun\n   core: %d x %d x %d\n '...
        'length: %-d, %d, %-d\n '...
        domainFormatString ' vertical scale = %-.2g\n'];
fprintf(str, r1, r2, r3, r1, r2, r3, m, n, p, ...
        dom(1), dom(2), dom(3), dom(4), dom(5), dom(6), vscl);

if ( loose )
    fprintf('\n');
end

end