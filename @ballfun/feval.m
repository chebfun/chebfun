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
            % Test 1:
            R_ndgrid = reshape(repmat(R(:,1,1),n*p,1),m,n,p);
            Lam_ndgrid = repmat(kron(reshape(Lam(1,:,1),n,1),ones(m,1)),p,1);
            Th_ndgrid = kron(reshape(Th(1,1,:),p,1),ones(m*n,1));
            ndgrid_test1 = ( norm( R(:) - R_ndgrid(:), inf) == 0 ) && ...
                ( norm( Lam(:) - Lam_ndgrid(:), inf) == 0 ) && ...
                ( norm( Th(:) - Th_ndgrid(:), inf) == 0 );
            
            % Test 2:
            R_ndgrid = reshape(repmat(R(:,1,1),n*p,1),m,n,p);
            Lam_ndgrid = kron(reshape(Lam(1,1,:),p,1),ones(m*n,1));
            Th_ndgrid = repmat(kron(reshape(Th(1,:,1),n,1),ones(m,1)),p,1);
            ndgrid_test2 = ( norm( R(:) - R_ndgrid(:), inf) == 0 ) && ...
                ( norm( Lam(:) - Lam_ndgrid(:), inf) == 0 ) && ...
                ( norm( Th(:) - Th_ndgrid(:), inf) == 0 );
            
            % Test 3:
            R_ndgrid = repmat(kron(reshape(R(1,:,1),n,1),ones(m,1)),p,1);
            Lam_ndgrid = kron(reshape(Lam(1,1,:),p,1),ones(m*n,1));
            Th_ndgrid = reshape(repmat(Th(:,1,1),n*p,1),m,n,p);
            ndgrid_test3 = ( norm( R(:) - R_ndgrid(:), inf) == 0 ) && ...
                ( norm( Lam(:) - Lam_ndgrid(:), inf) == 0 ) && ...
                ( norm( Th(:) - Th_ndgrid(:), inf) == 0 );
            
            % Test 4:
            R_ndgrid = repmat(kron(reshape(R(1,:,1),n,1),ones(m,1)),p,1);
            Th_ndgrid = kron(reshape(Th(1,1,:),p,1),ones(m*n,1));
            Lam_ndgrid = reshape(repmat(Lam(:,1,1),n*p,1),m,n,p);
            ndgrid_test4 = ( norm( R(:) - R_ndgrid(:), inf) == 0 ) && ...
                ( norm( Lam(:) - Lam_ndgrid(:), inf) == 0 ) && ...
                ( norm( Th(:) - Th_ndgrid(:), inf) == 0 );
            
            % Test 5:
            R_ndgrid = kron(reshape(R(1,1,:),p,1),ones(m*n,1));
            Th_ndgrid = repmat(kron(reshape(Th(1,:,1),n,1),ones(m,1)),p,1);
            Lam_ndgrid = reshape(repmat(Lam(:,1,1),n*p,1),m,n,p);
            ndgrid_test5 = ( norm( R(:) - R_ndgrid(:), inf) == 0 ) && ...
                ( norm( Lam(:) - Lam_ndgrid(:), inf) == 0 ) && ...
                ( norm( Th(:) - Th_ndgrid(:), inf) == 0 );
            
            % Test 6:
            R_ndgrid = kron(reshape(R(1,1,:),p,1),ones(m*n,1));
            Lam_ndgrid = repmat(kron(reshape(Lam(1,:,1),n,1),ones(m,1)),p,1);
            Th_ndgrid = reshape(repmat(Th(:,1,1),n*p,1),m,n,p);
            ndgrid_test6 = ( norm( R(:) - R_ndgrid(:), inf) == 0 ) && ...
                ( norm( Lam(:) - Lam_ndgrid(:), inf) == 0 ) && ...
                ( norm( Th(:) - Th_ndgrid(:), inf) == 0 );
            
            % Relate to fevalm code for efficiency:
            if ( ndgrid_test1 )
                
                % [R, Lam, Th] = ndgrid( r, lam, th):
                out = fevalm( f, R(:,1,1), reshape(Lam(1,:,1),n,1), reshape(Th(1,1,:),p,1) );
                
            elseif ( ndgrid_test2 )
                
                % [R, Th, Lam] = ndgrid( r, th, lam): 
                out = fevalm(f, R(:,1,1), reshape(Lam(1,1,:),p,1), reshape(Th(1,:,1),n,1));
                out = permute( out, [1 3 2]);
                
            elseif ( ndgrid_test3 )
                
                % [Th, R, Lam] = ndgrid( th, r, lam):
                out = fevalm(f, reshape(R(1,:,1),n,1), reshape(Lam(1,1,:),p,1), Th(:,1,1));
                out = permute( out, [3 1 2]);
                
            elseif ( ndgrid_test4 )
                
                % [Lam, R, Th] = ndgrid( lam, r, th);
                out = fevalm(f, reshape(R(1,:,1),n,1), Lam(:,1,1), reshape(Th(1,1,:),p,1));
                out = permute( out, [2 1 3]);
                
            elseif ( ndgrid_test5 )
                
                % [Lam, Th, R] = ndgrid( lam, th, r);
                out = fevalm(f, reshape(R(1,1,:),p,1), Lam(:,1,1), reshape(Th(1,:,1),n,1));
                out = permute( out, [2 3 1]);
                
            elseif ( ndgrid_test6 )
                
                % [Th, Lam, R] = ndgrid( th, lam, r);
                out = fevalm(f, reshape(R(1,1,:),p,1), reshape(Lam(1,:,1),n,1), Th(:,1,1));
                out = permute( out, [3 2 1]);
                
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
        
    else
        error('CHEBFUN:BALLFUN:feval:inputdim', ...
            'Evaluation points should all be the same dimension');
    end
else
    error('CHEBFUN:BALLFUN:feval:inputs', ...
        'Unrecognized arguments for evaluation.');
end
end