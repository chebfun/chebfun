function Y = harmonic(L, m, type)
%HARMONIC   Normalized, real-valued, eigenfunction of the Laplacian 
%           on the disk.
%
%   Y = HARMONIC(L, M) returns the eigenfunction  V_L^M(t,r) for
%   homogeneous Dirichlet boundary conditions. Here,
%        -pi <= t <= pi  is the angular coordinate, and
%          0 <= r  <= 1  is the radial coordinate.
%
%   Y = HARMONIC(L, M, 'dirichlet') is the same as Y = HARMONIC(L, M).
%
%   Y = HARMONIC(L, M,'neumann') returns the eigenfunction V_L^M(t,r)
%   for homogeneous Neumann boundary conditions.
%
%   The harmonic functions are the eigenfunctions of the  Laplacian on 
%   the disk, and they form an orthogonal basis with respect to the polar 
%   measure. For Dirichlet boundary conditions and L selected >= 0, they 
%   are of the form V_L^M(t, r) = A*cos(iLt).*J_L(a_Mr), where J_L is the Lth
%   j-Bessel function, a_M is the Mth positive zero of the function, and A
%   is a normalization factor. When L is specified as a negative value, 
%   they are of the form V_L^M(t, r) = A*sin(i|L|t).*J_|L|(a_Mr). 
%
% See also SPHEREFUN/SPHHARM.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 3 || isempty(type) )
    type = 'dirichlet';
end

Lsign = sign( L+1 );
L = abs( L );


if ( strcmpi(type, 'dirichlet') )
    % Calculate the mth positive zero of the Lth Bessel function
    % We would like to do the following (which is much faster). However,
    % that algorithm is not particularly accurate. 
    %  Jzero = besselroots(L, m);
    %  Jzero = Jzero( m );
    
    % Instead, we currently use chebfun/roots:
    Jzero = roots(chebfun(@(x) besselj(L,x), [sqrt((3/4)^2*pi^2+L^2) (m+L/2)*pi]));
    Jzero = Jzero(m);
    
    % Compute the normalization
    nrmlz = sqrt(2)/(sqrt((1+double(L==0))*pi)*abs(besselj(L+1, Jzero)));

    % Create a diskfun object for cylindrical harmonic:
    if ( abs( max(0, Lsign) ) == 1 )
        % Use cosine(L*theta) basis
        Y = diskfun(@(theta, r) cos(L*theta).*besselj(L, r*Jzero)*nrmlz, 'polar');
    else
        % Use sin(L*theta) basis
        Y = diskfun(@(theta, r) sin(L*theta).*besselj(L, r*Jzero)*nrmlz, 'polar');
    end
elseif ( strcmpi(type, 'neumann') )    
    if ( L==0 && m==1 )
        % Case where eigenvalue is zero; constant mode
        Y = diskfun(1);
    else
        % The x-term is not present when L=0
        Jzero = roots(chebfun(@(x) L*besselj(L,x)-x.^(min(1,L))...
            .*besselj(L+1, x), [L (m+L/2)*pi])); 
        Jzero = Jzero(m); % When L=0 this is still the correct index
        
        % Compute the normalization
        nrmlz = sqrt(2)/(sqrt((1-L^2/Jzero^2)*(1+double(L==0))*pi)*abs(besselj(L,Jzero)));

        % Create a diskfun object for the harmonic:
        if ( abs( max(0, Lsign) ) == 1 )
            % Use cosine(L*theta) basis
            Y = diskfun(@(theta, r) cos(L*theta).*besselj(L, r*Jzero)*nrmlz, 'polar');
        else
            % Use sin(L*theta) basis
            Y = diskfun(@(theta, r) sin(L*theta).*besselj(L, r*Jzero)*nrmlz, 'polar');
        end
    end   
else    
    error('CHEBFUN:DISKFUN:HARMONIC:TYPE', ['Unrecognized harmonic type.'...
        ' Only dirichlet and neumann are available types.'] )
end

end

% The old way of constructing cylindrical harmonics was to wrap this code
% in an annonymous function, but that was slow.  We are keeping this code
% here for historical purposes.
% function Z = cylindricalharmonic(theta, r, L, m, type)
% %CYLINDRICALHARMONIC   Construct function handle for cylindrical harmonic.
% %
% %   Z = CYLINDRICALHARMONIC( THETA, R, L, M, TYPE ), construct a function
% %   handle in coordinates (THETA, R) for the (L,M) cylindrical harmonic of
% %   type TYPE. If TYPE = 'dirichlet', then the harmonic has homogeneous 
% %   Dirichlet data; otherwise, TYPE = 'neumann', then the harmonic has 
% %   homogeneous neumann data at the boundary. 
% 
% % Calculate the mth positive zero of the Lth Bessel function
% Lsign = sign( L+1 );
% L = abs( L );
% 
% if ( strcmpi(type, 'dirichlet') )
%     % We would like to do the following (which is much faster). However,
%     % that algorithm is not particularly accurate. 
%     %  Jzero = besselroots(L, m);
%     %  Jzero = Jzero( m );
%     
%     % Instead, we currently use chebfun/roots:
%     Jzero = roots(chebfun(@(x) besselj(L,x), [sqrt((3/4)^2*pi^2+L^2) (m+L/2)*pi]));
%     Jzero = Jzero(m);
%     z = besselj(L, r*Jzero);
% 
%     % Choose sin or cos based on L value
%     pos = abs( max(0, Lsign) );
%     Z = ((pos)*cos(L*theta) + (1-pos)*sin(L*theta)).*z;
%     
%     % Normalize
%     k = double( L==0 );
%     Z = sqrt(2)/ (sqrt((1+k)*pi)*abs(besselj(L+1, Jzero)))*Z;
% elseif ( strcmpi(type, 'neumann') )    
%     if ( L==0 && m==1 )
%         % Case where eigenvalue is zero; constant mode
%         Z = 1 + r*0;  
%     else
%         % The x-term is not present when L=0
%         Jzero = roots(chebfun(@(x) L*besselj(L,x)-x.^(min(1,L))...
%             .*besselj(L+1, x), [L (m+L/2)*pi])); 
%         Jzero = Jzero(m); % When L=0 this is still the correct index
%         z = besselj(L, r*Jzero);
%         % Choose sin or cos based on L value
%         pos = abs(max(0,Lsign));
%         Z = ((pos)*cos(L*theta) +(1-pos)*sin(L*theta)).*z;
%         % Normalize
%         k = double(L==0);
%         Z = sqrt(2)/ (sqrt( (1-L^2/Jzero^2)*(1+k)*pi)*abs(besselj(L, Jzero)))*Z;
%     end   
% else    
%     error('DISKFUN:HARMONIC:TYPE', 'Unrecognized Bessel type.')
% end
% 
% end