function out = roots(f, varargin)
%ROOTS	Roots of a FUNCHEB in the interval [-1,1].
%   ROOTS(F) returns the real roots of the FUNCHEB F in the interval [-1,1].
%
%   ROOTS(F, PROP1, VAL1, PROP2, VAL2, ...) modifies the default ROOTS
%   properties. The PROPs and VALs may be any of the following:
%
%   ALL: 
%       [0] - Return only real-valued roots in [-1,1]
%        1  - Return roots outside of [-1,1] (including complex roots).
%
%   RECURSE:
%        0  - Compute roots without bisection (slower).
%       [1] - Bisect until length(F) < 50. (fast, but additional complex roots).
%
%   PRUNE:
%       [0]
%        1  - Prune 'spurious' complex roots if ALL == 1 and RECURSE == 0.
%
%   HSCALE:
%       [1] - Horizontal scale for adjusting relative tolerances.
%     double
%
%   If F is a vector-valued FUNCHEB then there is no reason to expect each
%   column to have the same number of roots. In order to return a useful output,
%   the roots of each column are computed and then padded with NaNs so that a
%   matrix may be returned. The columns of R = ROOTS(F) correspond to the
%   columns of F.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROOTS works by recursively bisecting the interval until the resulting funcheb
% is of degree less than 50, at which point a companion matrix is constructed to
% compute the roots.
%
% ROOTS performs all operations in coefficient space. In this representation,
% two matrices, Tleft and Tright, are constructed such that Tleft*c and Tright*c
% are the coefficients of the polynomials in the left and right intervals
% respectively. This is faster than evaluating the polynomial using barycentric
% interpolation in the respective intervals despite both computations requiring
% O(N^2) operations.
%
% For polynomials of degree larger than 512, the interval is bisected by
% evaluating on the left and right intervals using the Clenshaw algorithm.
%
% [Mathematical references]:
%  * I. J. Good, "The colleague matrix, a Chebyshev analogue of the companion
%    matrix", Quarterly Journal of Mathematics 12 (1961).
%
%  * J. A. Boyd, "Computing zeros on a real interval through Chebyshev expansion
%    and polynomial rootfinding", SIAM Journal on Numerical Analysis 40 (2002).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Deal with empty case:
if ( isempty(f) )
    out = [];
    return
end

% Support for vectorised FUNCHEB objects.
if ( size(f.values, 2) > 1 )
    % Copy f into g:
    g = f;
    % Initialise a cell array to hold roots of each column:
    r = cell(1, size(f.values, 2));
    % Loop over the columns:
    for j = 1:size(f.values, 2)
        g.values = f.values(:,j);
        g.coeffs = f.coeffs(:,j);
        g.vscale = f.vscale(j);
        r{j} = roots(g, varargin{:}); 
    end
    % Find the max length of r:
    mlr = max(cellfun(@length, r)); 
    % Pad the columns in r with NaNs:
    r = cellfun(@(x) [x ; NaN(mlr-length(x), 1)], r, 'UniformOutput', false);
    % Convert to an array for output:
    out = cell2mat(r);
    return
end

% Default preferences:
rootspref = struct('all', 0, 'recurse', 1, 'prune', 0, 'hscale', 1);
splitpoint = -0.004849834917525;

% Filter out the arguments:
j = 1;
while ( j <= length(varargin) )
    if any(strcmp(lower(varargin{j}), fieldnames(rootspref))) %#ok<STCI>
        rootspref.(varargin{j}) = varargin{j+1};
        j = j + 2;
    elseif strcmpi(varargin{j}, 'complex')
        rootspref.all = varargin{j+1};
        j = j + 2;
    else
        j = j + 1;
    end
end

% Trivial case for f constant:
if ( length(f) == 1 )
    if ( f.values(1) == 0 ) % return a root at centre of domain
        out = 0;
    else
        out = [];
    end
    return
end

% Get scaled coefficients for the recursive call:
c = flipud(f.coeffs) / f.vscale;

hscale = rootspref.hscale;

% Call the recursive rootsunit function:
r = rootsunit_coeffs( c , 100*eps*max(hscale, 1) );

% Prune the roots, if required:
if ( rootspref.prune && ~rootspref.recurse )
    rho = sqrt(eps)^(-1/length(f));
    rho_roots = abs(r + sqrt(r.^2 - 1));
    rho_roots(rho_roots < 1) = 1./rho_roots(rho_roots < 1);
    out = r(rho_roots <= rho);
else
    out = r;
end

    function r = rootsunit_coeffs ( c , htol  )
    % Computes the roots of the polynomial given by the coefficients
    % c on the unit interval.

        % Define these as persistent, need to compute only once.
        persistent Tleft Tright;

        % Simplify the coefficients
        n = length(c);
        tailMmax = 1e-15*norm(c,1);
        while ( (n > 1) && (abs(c(n)) <= tailMmax) )
            n = n - 1; 
        end
        
        % Wrap (i.e., alias), don't just truncate.
        if ( n > 1 && n < length(c) )
            c = funcheb.alias(c(end:-1:1), n);
            c = c(end:-1:1);
        end

        % Trivial case, n == 1
        if ( n == 1 )

            % If the function is zero, then place a root in the middle
            if ( c(1) == 0 )
                r = 0.0;
            else
                r = [];
            end

        % Trivial case, n == 2
        elseif ( n == 2 )

            % is the root in [-1,1]?
            r = -c(1) / c(2);
            if ( ~rootspref.all )
                if ( abs(imag(r)) > htol ) || ( r < -(1+htol) ) || ( r > (1+htol) )
                    r = [];
                else
                    r = max( min( real(r) , 1 ) , -1 );
                end
            end

        % Is n small enough to compute the roots directly?
        elseif ( ~rootspref.recurse || ( n <= 50 ) )

            % adjust the coefficients for the colleague matrix
            c = -0.5 * c(1:end-1) / c(end);
            c(end-1) = c(end-1) + 0.5;
            oh = 0.5 * ones(length(c)-1,1);

            % Modified colleague matrix:
            A = diag(oh,1) + diag(oh,-1);
            A(end-1,end) = 1;
            A(:,1) = flipud(c);

            % compute roots as eig(A)
            r = eig(A);

            % Clean the roots up a bit
            if ( ~rootspref.all )
            
                % Remove dangling imaginary parts
                mask = abs(imag(r)) < htol;
                r = real( r(mask) );
                
                % keep roots inside [-1 1]
                r = sort( r(abs(r) <= 1+2*htol) );
                
                % Correct roots over ends
                if ~isempty(r)
                    r(1) = max(r(1),-1);
                    r(end) = min(r(end),1);
                end

            % Prune?
            elseif ( rootspref.prune )
                rho = sqrt(eps)^(-1/n);
                rho_roots = abs(r+sqrt(r.^2-1));
                rho_roots(rho_roots < 1) = 1./rho_roots(rho_roots<1);
                r = r(rho_roots <= rho);
            end
            
        % Can we compute the new coefficients with a cheap matrix-vector?
        elseif ( n <= 513 )

            % Have we assembled the matrices Tleft and Tright?
            if isempty( Tleft )

                % create the coefficients for Tleft using the fft directly.
                x = chebptsAB(513, [-1, splitpoint]);
                Tleft = ones(513); 
                Tleft(:,2) = x;
                for k = 3:513
                    Tleft(:,k) = 2 * x .* Tleft(:,k-1) - Tleft(:,k-2); 
                end
                Tleft = [ Tleft(513:-1:2,:) ; Tleft(1:512,:) ];
                Tleft = real( fft( Tleft ) / 512 );
                Tleft = triu( [ 0.5*Tleft(1,:) ; Tleft(2:512,:) ; 0.5*Tleft(513,:) ] );

                % create the coefficients for Tright much in the same way
                x = chebptsAB(513, [splitpoint,1]);
                Tright = ones(513); 
                Tright(:,2) = x;
                for k = 3:513
                    Tright(:,k) = 2 * x .* Tright(:,k-1) - Tright(:,k-2); 
                end
                Tright = [ Tright(513:-1:2,:) ; Tright(1:512,:) ];
                Tright = real( fft( Tright ) / 512 );
                Tright = triu( [ 0.5*Tright(1,:) ; Tright(2:512,:) ; 0.5*Tright(513,:) ] );

            end % isempty(Tleft)

            % compute the new coefficients
            cleft = Tleft(1:n,1:n) * c;
            cright = Tright(1:n,1:n) * c;

            % recurse
            r = [ (splitpoint-1)/2 + (splitpoint+1)/2*rootsunit_coeffs( cleft , 2*htol )
                  (splitpoint+1)/2 + (1-splitpoint)/2*rootsunit_coeffs( cright , 2*htol ) ];

        % Otherwise, split using more traditional methods
        else
            
            % Evaluate the polynomial on both intervals
            v = funcheb.clenshaw( [ chebptsAB(n, [-1, splitpoint]) ; chebptsAB(n, [splitpoint, 1]) ], c(end:-1:1) );

            % Get the coefficients on the left
            cleft = f.chebpoly(v(1:n));
            cleft = cleft(end:-1:1);

            % Get the coefficients on the right
            cright = f.chebpoly(v(n+1:end));
            cright = cright(end:-1:1);

            % recurse
            r = [ (splitpoint-1)/2 + (splitpoint+1)/2*rootsunit_coeffs( cleft , 2*htol )
                  (splitpoint+1)/2 + (1-splitpoint)/2*rootsunit_coeffs( cright , 2*htol ) ];

        end

    end

    function y = chebptsAB(n, ab)
    % Y = CHEBPTSAB(N, [A, B]) is the N-point Chebyshev grid mapped to [A,B].

        a = ab(1);
        b = ab(2);
        x = f.chebpts(n);                 % [-1,1] grid
        y = b*(x+1)/2 + a*(1-x)/2;        % new grid

    end

end



