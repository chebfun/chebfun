function out = feval(varargin)
%FEVAL   Evaluate a BALLFUN at one or more points.
%   Y = FEVAL( F, X, Y, Z) evaluates a BALLFUN F at a point (X,Y,Z) in Cartesian
%   coordinates, where X, Y and Z are doubles.
%
%   Y = FEVAL( F, R, LAM, TH, 'spherical') evaluates a BALLFUN F in
%   spherical coordinates (R,LAM,TH). Here R, LAM and THETA are doubles representing 
%   the radius, azimuthal and polar angles (in radians) and must be points 
%   in the unit ball.
%
% See also FEVALM. 

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.\

f = varargin{1};

% Empty check:
if ( isempty(f) )
    out = [];
    return
end

% Figure out if Cartesian or polar coordinates should be used.
% Search for user-supplied 'spherical' flag in arguments:
isPolar = any(find(strcmp(varargin,'polar'))) || any(find(strcmp(varargin,'spherical')));
if ( ~isPolar )
X = varargin{2};
Y = varargin{3};
Z = varargin{4};
[Lam,Th,R] = cart2sph(X,Y,Z);
Th = pi/2 - Th;
if ( any(R > 1+1e-8) ) % Check for points off ball
    error('CHEBFUN:BALLFUN:FEVAL:pointsNotOnDisk',...
        ['The specified points to evaluate the function do not '...
        'lie sufficiently close to the unit ball.']);
end
out = feval(f, R, Lam, Th, 'spherical');

else
R = varargin{2};
Lam = varargin{3};
Th = varargin{4};
if ( isnumeric(R) && isnumeric( Lam ) && isnumeric( Th ) )
    if ( isequal(size(R), size(Lam)) && isequal(size(R), size(Th)) )
        if ( isscalar(R) && isscalar(Lam) && isscalar(Th) )
            out = fevalm( f, R, Lam, Th);
        else
            
            % If the evaluation points are derived from ndgrid, then there is a
            % fast way to evaluate a BALLFUN. Check for this property.
            [m, n, p] = size(R);
            % This is the most common ndgrid input to BALLFUN, but not the only
            % one.
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
                            out(i1, j1, k1) = feval(f, R(i1,j1,k1), Lam(i1,j1,k1), Th(i1,j1,k1), 'spherical');
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

% Return real values if the function is real
if f.isReal
    out = real(out);
end
end
