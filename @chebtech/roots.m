function out = roots(f, varargin)
%ROOTS   Roots of a CHEBTECH in the interval [-1,1].
%   ROOTS(F) returns the real roots of the CHEBTECH F in the interval [-1,1].
%
%   ROOTS(F, PROP1, VAL1, PROP2, VAL2, ...) modifies the default ROOTS
%   properties. The PROPs (strings) and VALs may be any of the following:
%
%   ALL: 
%       [0] - Return only real-valued roots in [-1,1].
%        1  - Return roots outside of [-1,1] (including complex roots).
%
%   COMPLEX:
%       [0] - No effect.
%        1  - Equivalent to setting both PRUNE and ALL = 1.
%
%   FILTER:
%       [ ]
%   @filter(R,F) - A function handle which accepts the sorted computed roots, R, 
%                  and the CHEBTECH, F, and filters the roots as it see fit.
%   RECURSE:
%        0  - Compute roots without interval subdivision (slower).
%       [1] - Subdivide until length(F) < 50. (Can cause additional complex
%             roots).
%
%   PRUNE:
%       [0]
%        1  - Prune 'spurious' complex roots if ALL == 1 and RECURSE == 0.
%
%   QZ: 
%       [0] - Use the colleague matrix linearization and the QR algorithm.
%        1  - Use the colleague matrix pencil linearization and the QZ 
%             algorithm for potentially extra numerical stability. 
%
%   ZEROFUN:
%        0  - Return empty if F is identically 0.
%       [1] - Return a root at x = 0 if F is identically 0.
%
%   If F is an array-valued CHEBTECH then there is no reason to expect each
%   column to have the same number of roots. In order to return a useful output,
%   the roots of each column are computed and then padded with NaNs so that a
%   matrix may be returned. The columns of R = ROOTS(F) correspond to the
%   columns of F.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROOTS works by recursively subdividing the interval until the resulting
% CHEBTECH is of degree less than 50, at which point a colleague matrix is
% constructed to compute the roots.
%
% ROOTS performs all operations in coefficient space. In this representation,
% two matrices, TLEFT and TRIGHT (both of size 512 by 512), are constructed such
% that TLEFT*C and TRIGHT*C are the coefficients of the polynomials in the left
% and right intervals respectively. This is faster than evaluating the
% polynomial using the barycentric formula or Clenshaw's algorithm in the
% respective intervals despite both computations requiring O(N^2) operations.
%
% For polynomials of degree larger than 512, the interval is subdivided by
% evaluating on the left and right intervals using the Clenshaw algorithm. The
% subdivision occurs at an arbitrary point _near_ but not _at_ the centre of the
% domain (in fact, -0.004849834917525 on [-1 1]) to avoid introducing additional
% spurious roots (since x = 0 is often a special point).
%
% Note that ROOTS uses CHEBTECH2 technology to subdivide the interval,
% regardless of whether F is a CHEBTECH1 or a CHEBTECH2.
%
% [Mathematical references]:
%  * I. J. Good, "The colleague matrix, a Chebyshev analogue of the companion
%    matrix", Quarterly Journal of Mathematics 12 (1961).
%
%  * J. A. Boyd, "Computing zeros on a real interval through Chebyshev expansion
%    and polynomial rootfinding", SIAM Journal on Numerical Analysis 40 (2002).
%
%  * L. N. Trefethen, Approximation Theory and Approximation Practice, SIAM,
%    2013, chapter 18.
%
%  [TODO]: Update this reference.
%  * Y. Nakatsukasa and V. Noferini, On the stability of polynomial rootfinding
%  via linearizations in nonmonomial bases, (In Prep).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Deal with empty case:
if ( isempty(f) )
    out = [];
    return
end

if ( size(f.coeffs, 2) == 1 )   % F is a scalar-value CHEBTECH.
    
    % Simply call roots_scalar():
    out = roots_scalar(f, varargin{:});
    
else                            % Support for array-valued CHEBTECH objects.

    % Initialise a cell array to hold roots of each column:
    r = cell(1, size(f.coeffs, 2));
    
    % Copy columns of f into an array g:
    g = mat2cell(f);

    % Loop over the columns of f / elements of g:
    for j = 1:size(f.coeffs, 2)
        r{j} = roots_scalar(g{j}, varargin{:}); 
    end

    % Find the max length of r:
    mlr = max(cellfun(@length, r)); 

    % Pad the columns in r with NaNs:
    r = cellfun(@(x) [x ; NaN(mlr - length(x), 1)], r, 'UniformOutput', false);

    % Convert to an array for output:
    out = cell2mat(r);

end    

end

function out = roots_scalar(f, varargin)

% Default preferences:
rootsPref = struct('all', 0, 'recurse', 1, 'prune', 0, 'zeroFun', 1, ...
    'qz', 0, 'filter', []);
% Subdivision maps [-1,1] into [-1, splitPoint] and [splitPoint, 1].
splitPoint = -0.004849834917525; % This is an arbitrary number.

if ( nargin > 1 && isa(varargin{1}, 'struct') )
    rootsPref = varargin{1};
    varargin(1) = [];
end

% Filter out the arguments:
j = 1;
while ( j < length(varargin) )
    if ( strcmpi(varargin{j}, 'complex') && varargin{j+1} )
        rootsPref.prune = true;
        rootsPref.all = true;
        j = j + 2;
    elseif ( strcmpi(varargin{j}, 'all') )
        rootsPref.all = varargin{j+1};
    elseif ( strcmpi(varargin{j}, 'recurse') )
        rootsPref.recurse = varargin{j+1};
        j = j + 2;    
    elseif ( strcmpi(varargin{j}, 'prune') )
        rootsPref.prune = varargin{j+1};
        j = j + 2;            
    elseif ( strcmpi(varargin{j}, 'zeroFun') )
        rootsPref.zeroFun = varargin{j+1};        
        j = j + 2;
    elseif ( strcmpi(varargin{j}, 'qz') )
        rootsPref.qz = varargin{j+1};
        j = j + 2;        
    elseif ( strcmpi(varargin{j}, 'filter') )
        rootsPref.filter = varargin{j+1};
        j = j + 2;                
    else
        j = j + 1;
    end
end

% Trivial case for f constant:
if ( length(f) == 1 )
    if ( f.coeffs(1) == 0 && rootsPref.zeroFun )
        % Return a root at centre of domain:
        out = 0;
    else
        % Return empty:
        out = [];
    end
    return
end

% Get scaled coefficients for the recursive call:
c = f.coeffs/vscale(f);

% Call the recursive rootsunit function:
% TODO:  Does the tolerance need to depend on some notion of hscale?
r = rootsunit_coeffs(c, 100*eps);

% Try to filter out spurious roots:
if ( ~isempty(rootsPref.filter) )
    r = sort(r, 'ascend');
    r = rootsPref.filter(r, f);
end

% Prune the roots, if required:
if ( rootsPref.prune && ~rootsPref.recurse )
    rho = sqrt(eps)^(-1/length(f));
    rho_roots = abs(r + sqrt(r.^2 - 1));
    rho_roots(rho_roots < 1) = 1./rho_roots(rho_roots < 1);
    out = r(rho_roots <= rho);
else
    out = r;
end

    function r = rootsunit_coeffs(c, htol)
    % Computes the roots of the polynomial given by the coefficients c on the
    % unit interval.

        % Define these as persistent, need to compute only once.
        persistent Tleft Tright

        % Simplify the coefficients:
        tailMmax = eps*norm(c, 1);
        % Find the final coefficient about tailMax:
        n = find(abs(c) > tailMmax, 1, 'last');

        % Should we alias or truncate here? We truncate here for speed (about
        % 30-50% faster on example with a large amount of subdivision. 
        % Wrap (i.e., alias), don't just truncate:
%         if ( ~isempty(n) && (n > 1) && (n < length(c)) )
%             c = chebtech2.alias(c(end:-1:1), n);
%             c = c(end:-1:1);
%         end
        % Truncate the coefficients (rather than alias):
        if ( ~isempty(n) && (n > 1) && (n < length(c)) )
            c = c(1:n);
        end

        % Trivial case, n == []:
        if ( isempty(n) )
            
            if ( rootsPref.zeroFun )
                % If the function is zero, then place a root in the middle:
                r = 0;
            else
                % Else return empty:
                r = [];
            end
            
        % Trivial case, n == 1:
        elseif ( n == 1 )

            % If the function is zero, then place a root in the middle:
            if ( c(1) == 0 && rootsPref.zeroFun )
                r = 0;
            else
                % Else return empty:
                r = [];
            end

        % Trivial case, n == 2:
        elseif ( n == 2 )

            % Is the root in [-1,1]?
            r = -c(1)/c(2);
            if ( ~rootsPref.all )
                if ( (abs(imag(r)) > htol) || ...
                     (r < -(1 + htol)) || ...
                     (r > (1 + htol)) )
                    r = [];
                else
                    r = max(min(real(r), 1), -1);
                end
            end

        % Is n small enough for the roots to be calculated directly?
        elseif ( ~rootsPref.recurse || (n <= 50) )

            % Adjust the coefficients for the colleague matrix:
            cOld = c(:); 
            c = -0.5 * c(1:end-1) / c(end);
            c(end-1) = c(end-1) + 0.5;
            
            % Modified colleague matrix:
            % [TODO]: Would the upper-Hessenberg form be better?
            oh = 0.5 * ones(length(c)-1, 1);
            A = diag(oh, 1) + diag(oh, -1);
            A(end-1, end) = 1;
            A(:, 1) = flipud( c );
            
            % Compute roots by an EP if qz is 'off' and by a GEP if qz is 'on'.
            % QZ has been proved to be stable, while QR is not (see [Nakatsukasa
            % & Noferini, 2014]):
            if ( isfield(rootsPref, 'qz') && rootsPref.qz )
                % Set up the GEP. (This is more involved because we are scaling
                % for extra stability.)
                B = eye( size( A ) ); 
                cOld = cOld / norm( cOld, inf ); 
                B(1, 1) = cOld( end );
                cOld = -0.5 * cOld( 1:end-1 );
                cOld( end-1 ) = cOld( end-1 ) + 0.5 * B(1, 1);
                A(:, 1) = flipud( cOld ); 
                r = eig(A, B);
            else
                % Standard colleague (See [Good, 1961]):
                r = eig(A); 
            end
           
            % Clean the roots up a bit:
            if ( ~rootsPref.all )
                % Remove dangling imaginary parts:
                mask = abs(imag(r)) < htol;
                r = real( r(mask) );
                % Keep roots inside [-1 1]:
                r = sort( r(abs(r) <= 1 + htol) );
                % Correct roots over ends:
                if ( ~isempty(r) )
                    r(1) = max(r(1), -1);
                    r(end) = min(r(end), 1);
                end
            elseif ( rootsPref.prune )
                rho = sqrt(eps)^(-1/n);
                rho_roots = abs(r + sqrt(r.^2 - 1));
                rho_roots(rho_roots < 1) = 1./rho_roots(rho_roots < 1);
                r = r(rho_roots <= rho);
            end

        % If n <= 513 then we can compute the new coefficients with a
        % matrix-vector product.
        elseif ( n <= 513 )
            
            % Have we assembled the matrices TLEFT and TRIGHT?
            if ( isempty(Tleft) )
                % Create the coefficients for TLEFT using the FFT directly:
                x = chebptsAB(513, [-1, splitPoint]);
                Tleft = ones(513); 
                Tleft(:,2) = x;
                for k = 3:513
                    Tleft(:,k) = 2 * x .* Tleft(:,k-1) - Tleft(:,k-2); 
                end
                Tleft = [ Tleft(513:-1:2,:) ; Tleft(1:512,:) ];
                Tleft = real(fft(Tleft) / 512);
                Tleft = triu( [ 0.5*Tleft(1,:) ; Tleft(2:512,:) ; 0.5*Tleft(513,:) ] );

                % Create the coefficients for TRIGHT much in the same way:
                x = chebptsAB(513, [splitPoint,1]);
                Tright = ones(513); 
                Tright(:,2) = x;
                for k = 3:513
                    Tright(:,k) = 2 * x .* Tright(:,k-1) - Tright(:,k-2); 
                end
                Tright = [ Tright(513:-1:2,:) ; Tright(1:512,:) ];
                Tright = real(fft(Tright) / 512);
                Tright = triu( [ 0.5*Tright(1,:) ; Tright(2:512,:) ; 0.5*Tright(513,:) ] );
            end

            % Compute the new coefficients:
            cleft = Tleft(1:n,1:n) * c;
            cright = Tright(1:n,1:n) * c;

            % Recurse:
            r = [ (splitPoint - 1)/2 + (splitPoint + 1)/2*rootsunit_coeffs(cleft, 2*htol) ;
                  (splitPoint + 1)/2 + (1 - splitPoint)/2*rootsunit_coeffs(cright, 2*htol) ];

        % Otherwise, split using more traditional methods (i.e., Clenshaw):
        else
            
            % Evaluate the polynomial on both intervals:
            v = chebtech.clenshaw([ chebptsAB(n, [ -1, splitPoint ]) ; ...
                chebptsAB(n, [ splitPoint, 1 ]) ], c);

            % Get the coefficients on the left:
            cleft = chebtech2.vals2coeffs(v(1:n));            

            % Get the coefficients on the right:
            cright = chebtech2.vals2coeffs(v(n+1:end));           

            % Recurse:
            r = [ (splitPoint - 1)/2 + (splitPoint + 1)/2*rootsunit_coeffs(cleft, 2*htol) ;
                  (splitPoint + 1)/2 + (1 - splitPoint)/2*rootsunit_coeffs(cright, 2*htol) ];

        end

    end

end

function y = chebptsAB(n, ab)
% Y = CHEBPTSAB(N, [A, B]) is the N-point Chebyshev grid mapped to [A,B].

    a = ab(1);
    b = ab(2);
    x = chebtech2.chebpts(n);          % [-1,1] grid
    y = b*(x + 1)/2 + a*(1 - x)/2;     % new grid

end
