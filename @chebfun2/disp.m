function disp(F)
%DISP   Display a CHEBFUN2 to the command line.
% 
% See also DISPLAY.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

loose = strcmp( get(0, 'FormatSpacing'), 'loose' );

% Get display style and remove trivial empty CHEBFUN2 case. 
if ( isempty(F) )
    fprintf('    empty chebfun2\n')
    if ( loose )
        fprintf('\n');
    end
    return
end

% Get information that we want to display:
dom = F.domain;                           % Domain
len = length(F);                          % Numerical rank
[xx, yy] = meshgrid(dom(1:2), dom(3:4));   
vals = feval(F, xx, yy ).';               % Corner values
vals = vals(:);
vscl = vscale(F);                         % vertical scale

% Check underlying tech is a TRIGTECH: 
techCol = get(F.cols.funs{1}, 'tech');
techRow = get(F.rows.funs{1}, 'tech');

% Display the information: 
if ( isa(techCol(), 'trigtech') && isa(techRow(), 'trigtech') )
    disp('   chebfun2 object  (trig)')
elseif ( isa(techRow(), 'trigtech') )
    disp('   chebfun2 object  (trig in x)')
elseif ( isa(techCol(), 'trigtech') )
    disp('   chebfun2 object  (trig in y)')
else
    disp('   chebfun2 object')
end
fprintf('       domain                 rank       corner values\n');
if ( isreal(vals) )
    fprintf('[%4.2g,%4.2g] x [%4.2g,%4.2g]   %6i     [%4.2g %4.2g %4.2g %4.2g]\n', ...
        dom(1), dom(2), dom(3), dom(4), len, vals);
else
    fprintf('[%4.2g,%4.2g] x [%4.2g,%4.2g]   %6i     [  complex values  ]\n', ...
        dom(1), dom(2), dom(3) , dom(4), len);
end
fprintf('vertical scale = %3.2g \n', vscl)

if ( loose )
    fprintf('\n');
end

end