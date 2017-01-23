function display(F)
%DISP   Display a CHEBFUN3T object to the command line.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isequal(get(0, 'FormatSpacing'), 'compact') )
	disp([inputname(1), ' =']);
else
	disp(' ');
	disp([inputname(1), ' =']);
	disp(' ');
end

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
[m, n, p] = size(F.coeffs);               % Size of coeffs
vscl = F.vscale;                          % vertical scale

disp('   chebfun3t object ')

if ( all(floor(dom) == dom) ) 
    % Corners of the domain are all integers.
    domainFormatString = 'domain: [%-d, %d] x [%-d, %d] x [%-d, %d]';
    
else    
    domainFormatString = ['domain: [%-.3g, %.3g] x [%-.3g, %.3g] x '...
        '[%-.3g, %.3g]'];
end
    
str = ['   coeffs: %d x %d x %d \n   '...
    domainFormatString '\n   vertical scale = %-.2g\n'];
        
fprintf(str, m, n, p, dom(1), dom(2), dom(3), dom(4), dom(5), dom(6), vscl);

if ( loose )
    fprintf('\n');
end

end