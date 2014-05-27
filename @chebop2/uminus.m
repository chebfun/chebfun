function N = uminus(N)
%UMINUS unitary minus for chebop2 objects

A = N.coeffs; 
if ( iscell(A) )
    for jj = 1:size(A,1)
        for kk = 1:size(A,2) 
            A{jj,kk} = -A{jj,kk};
        end
    end
else
    A = -A; 
end


if ~isempty(N.lbc) || ~isempty(N.rbc) || ~isempty(N.ubc) || ~isempty(N.dbc)
    warning('Operator has BCs. These were not negated.');
end

N.coeffs = A; 
op = N.op; 
N.op = @(u) -op(u); 

end