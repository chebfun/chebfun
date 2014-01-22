function display(L)
%DISPLAY    Pretty-print a linop.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Obtain size of the linop
[m, n] = size(L);

% Determine whether format is compact or loose
loose = ~isequal(get(0,'FormatSpacing'),'compact');

if loose
   fprintf('\n'); 
end

% Print input variable name
disp([inputname(1) ' = ']);

if loose
   fprintf('\n'); 
end

fprintf('   %ix%i linear operator', m,n)

if loose
   fprintf('\n'); 
end

nc = size(L.constraint.operator,1);
if ( nc == 1 )
    fprintf('\n   with 1 constraint/boundary condition\n')
elseif ( nc > 0 )
    fprintf('\n   with %i constraints/boundary conditions\n',nc)
end

if loose
   fprintf('\n'); 
end

end
