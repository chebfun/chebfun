function disp(F)
%DISP   Display a CHEBFUN3 to the command line.
% 
%   See also DISPLAY.

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
[r1, r2, r3] = rank(F);                   % Rank
vscl = vscale(F);                         % vertical scale

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

if all(floor(dom) == dom) % Corners of the domain are all integers.
    fprintf(['   cols: [Inf x %d chebfun]\n   rows: [Inf x %d chebfun]\n  '...
        'tubes: [Inf x %d chebfun]\n   core: [%d x %d x %d double]\n '...
        'domain: [%-2d,%2d] x [%-2d,%2d] x [%-2d,%2d]\n vertical scale '...
        '= %-4.2g\n'], r1, r2, r3, r1, r2, r3, dom(1), dom(2), dom(3), ...
        dom(4), dom(5), dom(6), vscl);
else
    fprintf(['   cols: [Inf x %d chebfun]\n   rows: [Inf x %d chebfun]\n  '...
        'tubes: [Inf x %d chebfun]\n   core: [%d x %d x %d double]\n '...
        'domain: [%-3.3g, %3.3g] x [%-3.3g, %3.3g] x [%-3.3g, %3.3g]\n '...
        'vertical scale = %-4.2g\n'], r1, r2, r3, r1, r2, r3, dom(1), ...
        dom(2), dom(3), dom(4), dom(5), dom(6), vscl);
end

if ( loose )
    fprintf('\n');
end

end