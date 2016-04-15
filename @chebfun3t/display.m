function display(F)
%DISP   Display a CHEBFUN3T to the command line.

loose = strcmp( get(0, 'FormatSpacing'), 'loose');

% Get display style and remove trivial empty CHEBFUN3 case.
if ( isempty(F) )
    fprintf('    empty chebfun3t\n')
    if ( loose )
        fprintf('\n');
    end
    return
end

% Get information that we want to display:
dom = F.domain;                           % Domain
[m, n, p] = size(F.coeffs);             % Size of coeffs
vscl = F.vscale;                         % vertical scale

disp('   chebfun3t object ')

if all(floor(dom) == dom) % Corners of the domain are all integers.
fprintf(['   coeffs: [%d x %d x %d double]\n   '...
    'domain: [%-2d,%2d] x [%-2d,%2d] x [%-2d,%2d]\n vertical scale '...
    '= %-4.2g\n'], m, n, p, dom(1), dom(2), dom(3), ...
    dom(4), dom(5), dom(6), vscl);

else
    fprintf(['   coeffs: [%d x %d x %d double]\n   '...
    'domain: [%-3.3g, %3.3g] x [%-3.3g, %3.3g] x [%-3.3g, %3.3g]\n vertical scale '...
    '= %-4.2g\n'], m, n, p, dom(1), dom(2), dom(3), ...
    dom(4), dom(5), dom(6), vscl);
end

if ( loose )
    fprintf('\n');
end

end