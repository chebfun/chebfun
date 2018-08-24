function out = feval(f, R, Lam, Th)
%FEVAL Evaluate a BALLFUN.
%   FEVAL(F, R, LAM, TH) is the array of values of the BALLFUN
%   function F at the grid R x LAM x TH.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.\

% Empty check:
if ( isempty(f) )
    out = [];
    return
end

if ( isnumeric(R) && isnumeric( Lam ) && isnumeric( Th ) ) 
   if ( isequal(size(R), size(Lam)) && isequal(size(R), size(Th)) )
        if ( isscalar(R) && isscalar(Lam) && isscalar(Th) ) 
            out = fevalm( f, R, Lam, Th); 
        else
       
        % If the evaluation points are derived from ndgrid, then there is a
        % fast way to evaluate a BALLFUN. Check for this property.
        [m, n, p] = size(R);
        % This is the most common ndgrid input to BALLFUN, but not the only
        % one. TODO: Add the six different types of ndgrid inputs: 
        R_ndgrid = reshape(repmat(R(:,1,1),n*p,1),m,n,p);
        Lam_ndgrid = repmat(kron(reshape(Lam(1,:,1),n,1),ones(m,1)),p,1);
        Th_ndgrid = kron(reshape(Th(1,1,:),p,1),ones(m*n,1));
        ndgrid_test = ( norm( R(:) - R_ndgrid(:), inf) == 0 ) && ... 
                      ( norm( Lam(:) - Lam_ndgrid(:), inf) == 0 ) && ...
                      ( norm( Th(:) - Th_ndgrid(:), inf) == 0 );
        if ( ndgrid_test ) 
            % Just call fevalm code: 
            out = fevalm( f, R(:,1,1), reshape(Lam(1,:,1),n,1), reshape(Th(1,1,:),p,1) ); 
        else 
            % Evaluate at tensors, but didn't pass the ndgrid test.
            out = zeros(m, n, p);
            for i1 = 1:m 
                for j1 = 1:n 
                    for k1 = 1:p 
                        out(i1, j1, k1) = feval(f, R(i1,j1,k1), Lam(i1,j1,k1), Th(i1,j1,k1) ); 
                    end
                end
            end
        end
        end

   end
else
    error('CHEBFUN:BALLFUN:feval:inputs', ...
        'Unrecognized arguments for evaluation.');
end  
end