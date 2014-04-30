function display(F)
% DISPLAY   Display a CHEBFUN2V.
% 
% DISPLAY(F) outputs important information about the CHEBFUN2V F to the
% command window, including its domain of definition, length (number of 
% pivots used to represent it), and a summary of its structure. 
%
% It is called automatically when the semicolon is not used at the
% end of a statement that results in a CHEBFUN2V.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information. 

loose = strcmp(get(0,'FormatSpacing'),'loose');
if loose
    fprintf('\n%s = \n\n',inputname(1))
else
    fprintf('%s = \n',inputname(1))
end

% compact version
if ( isempty( F ) )
    fprintf('empty chebfun2v\n')
    return
end

if ( F.isTransposed )
    tString = 'Row vector';
else
    tString = 'Column vector';
end

disp(['chebfun2v object ' '(' tString ')' ])

% Display its two CHEBFUN2 halves.
for j = 1 : F.nComponents
    display( F.components{j} );
end

end
