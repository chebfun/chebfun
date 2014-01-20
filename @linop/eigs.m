function varargout = eigs(L,varargin)
% EIGS  Find selected eigenvalues and eigenfunctions of a linop.
% D = EIGS(A) returns a vector of 6 eigenvalues of the linop A. EIGS will
% attempt to return the eigenvalues corresponding to the least oscillatory
% eigenfunctions. (This is unlike the built-in EIGS, which returns the
% largest eigenvalues by default.)
%
% [V,D] = EIGS(A) returns a diagonal 6x6 matrix D of A's least oscillatory
% eigenvalues, and the corresponding eigenfunctions in V. If A operates on
% a single variable, then V is a quasimatrix of size Inf-by-6. If A
% operates on m (m>1) variables, then V is a 1-by-m cell array of
% quasimatrices. You can also use [V1,V2,...,Vm,D] = EIGS(A) to get
% a separate quasimatrix for each variable.
%
% EIGS(A,B) solves the generalized eigenproblem A*V = B*V*D, where B
% is another linop.
%
% EIGS(A,K) and EIGS(A,B,K) find the K smoothest eigenvalues.
%
% EIGS(A,K,SIGMA) and EIGS(A,B,K,SIGMA) find K eigenvalues. If SIGMA is a
% scalar, the eigenvalues found are the ones closest to SIGMA. Other
% selection possibilities for SIGMA are:
%    'LM' (or Inf) and 'SM' for largest and smallest magnitude
%    'LR' and 'SR' for largest and smallest real part
%    'LI' and 'SI' for largest and smallest imaginary part
% SIGMA must be chosen appropriately for the given operator. For
% example, 'LM' for an unbounded operator will fail to converge.
%
% Despite the syntax, this version of EIGS does not use iterative methods
% as in the built-in EIGS for sparse matrices. Instead, it uses the
% built-in EIG on dense matrices of increasing size, stopping when the
% targeted eigenfunctions appear to have converged, as determined by the
% chebfun constructor.
%
% EXAMPLE: Simple harmonic oscillator
%
%   d = domain(0,pi);
%   A = diff(d,2) & 'dirichlet';
%   [V,D] = eigs(A,10);
%   format long, sqrt(-diag(D))  % integers, to 14 digits
%
% See also EIGS, EIG.

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Parsing inputs.
M = [];  k = 6;  sigma = []; 
prefs = L.prefs;
discType = prefs.discretization;
gotk = false;
j = 1;
while (nargin > j)
    item = varargin{j};
    if isa(item,'linop')
        % Generalized operator term
        M = item;
%     elseif isstruct(varargin{j}) && isfield(varargin{j},'name')
%         if ~strcmp(varargin{j}.name,'linear')
%             map = varargin{j};
%         end
    elseif isa(item,'chebDiscretization')
        discType = item;
    else
        % k must be given before sigma.
        if ( ~gotk || ischar(item) )
            k = item;
            gotk = true;
        else
            sigma = item;
        end
    end
    j = j+1;
end

% Assign default to k if needed.
if ( isnan(k) || isempty(k) )
    k = 6; 
end

% maxdegree = cheboppref('maxdegree');
m = size(L,2);
if m ~= size(L,1)
    error('LINOP:eigs:notsquare','Block size must be square.')
end

%% Set up the discretisation:
if ( isa(discType, 'function_handle') )
    % Create a discretization object
    disc = discType(L);  
        
    % Set the allowed discretisation lengths: (TODO: A preference?)
    dimVals = L.prefs.dimensionValues;
    
    % Update the discretistion dimension on unhappy pieces:
    disc.dimension = repmat(dimVals(1), 1, numel(disc.domain)-1);
    dimVals(1) = [];
else
    % A discretisation is given:
    disc = discType;
    
    % Initialise dimVals;
    dimVals = max(disc.dimension);
end


dom = L.domain;
disc.domain = dom;

if ( isempty(L.continuity) )
     % Apply continuity conditions:
     disc.source = deriveContinuity(L);
end

discM = [];
if ( ~isempty(M) )
    dom = mergeDomains(disc,dom,M.domain);
    disc.domain = dom;
    dconstructor = str2fun( class(disc) );
    discM = dconstructor(M.disc.dimension,disc.domain);
    discM.source.constraint = disc.source.constraint;
    discM.source.continuity = disc.source.continuity;
end


numints = disc.numIntervals;
isFun = isFunVariable(L);

% 'SM' is equivalent to eigenvalues nearest zero.
if strcmpi(sigma,'SM')
    sigma = 0;
end

% Automatic mode: find the best sigma by going where the convergence appears to
% be fastest. 
if isempty(sigma)
    % Try to determine where the 'most interesting' eigenvalue is.
    disc.dimension = 33*ones(1,numints);
    [V1,D1] = getEigenvalues(disc,discM,33,0);
    disc.dimension(:) = 65;
    [V2,D2] = getEigenvalues(disc,discM,33,0);
    lam1 = diag(D1);  lam2 = diag(D2);
    dif = bsxfun(@minus,lam1.',lam2);
    delta = min( abs(dif) );   % diffs from 33->65
    bigdel = (delta > 1e-12*norm(lam1,Inf));
    
    % Trim off things that are still changing a lot (relative to new size).
    lam1b = lam1; lam1b(bigdel) = 0;
    bigdel = logical((delta > 1e-3*norm(lam1b,Inf)) + bigdel);
    
    if all(bigdel)
        % All values changed somewhat--choose the one changing the least.
        [~,idx] = min(delta);
        sigma = lam1(idx);
    else
        % One by one, convert the eigenvectors to functions and check their cheb
        % expansion coefficients.
        U = partition(disc,V2);  % each cell is array valued, for one variable
        
        % Combine the different variable components into a single variable for
        % coefficient conversion.
        Z = 0;
        for j = ( find(isFun) )
            Z = Z + U{j};
        end
        
        z = toFunction(disc,Z);
        
        % Compute the 1-norm of the polynomial expansions, summing over smooth
        % pieces, for all columns.
        onenorm = 0;
        for j = 1:disc.numIntervals
            onenorm = onenorm + sum( abs( chebpoly(z,j) ), 2 );
        end
        
        [~,index] = min(onenorm);  
        sigma = lam2(index);
    end
end



% Default settings
%settings = chebopdefaults;
%settings.scale = L.scale;

% Linear combination coefficients for convergence test.
coeff = 1./(2*(1:k)');

for dim = dimVals

    [V,D] = getEigenvalues(disc,discM,k,sigma);
        
    % Combine the eigenfunctions into a composite.
    v = V*coeff(1:size(V,2));
    
    % Convert the different components into cells
    u = partition(disc,v);
   
    % Test the happieness of the function pieces:
    [isDone, epsLevel] = testConvergence(disc, u(isFun));
    
    if ( all(isDone) )
        break
    else
        % Update the discretistion dimension on unhappy pieces:
        disc.dimension(~isDone) = dim;
    end
    
end

% Detect finite rank operators.
if size(D,1) < k
    if gotk
        warning('CHEBFUN:linop:eigs:rank',...
            'Input has finite rank, only %d eigenvalues returned.', size(D,1));
    end
    k = size(D,1);
end

if nargout < 2  % Return the eigenvalues only
    varargout = { diag(D) };
else            % Unwrap the eigenvectors for output  
    
    u = mat2fun(disc,V);
        
    % Find the norm in each eigenfunction (aggregated over variables).
    nrmsq = zeros(1,k);
    for j = 1:length(u)
        nrmsq = nrmsq + sum( u{j}.*conj(u{j}), 1 );
    end
    
    % Normalize each eigenfunction.
    scale = diag( 1./sqrt(nrmsq') );
    for j = 1:length(u)
        u{j} = u{j}*scale;
    end
    
     varargout = { chebmatrix(u), D };
end

end
% END OF MAIN FUNCTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [V,D] = getEigenvalues(disc,discM,k,sigma)
% Formulate the discrete problem and solve for the eigenvalues

    % Discretize the operator (incl. constraints/continuity):
    [A,P,C] = matrix(disc);
    nc = size(C,1);
    if ~isempty(discM)
        discM.dimension = disc.dimension;
        discM.domain = disc.domain;
        B = matrix(discM);
        B(1:nc,:) = 0;  % don't need the constraints on this side
    else
        B = [ zeros(nc,size(A,2)); P ];
    end
    
    if ( length(A) <= 2000 )
        [V,D] = eig(full(A),full(B));
        % Find the ones we're looking for.
        N = disc.dimension;
        idx = nearest(diag(D),V,sigma,min(k,N),N,disc);
        V = V(:,idx);
        D = D(idx,idx);
    else
        % TODO: Experimental.
        [V,D] = eigs(A,B,k,sigma);
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns index vector that sorts eigenvalues by the given criterion.
function idx = nearest(lam,V,sigma,k,N,disc)

if isnumeric(sigma)
    if isinf(sigma)
        [junk,idx] = sort(abs(lam),'descend');
    else
        [junk,idx] = sort(abs(lam-sigma));
    end
else
    switch upper(sigma)
        case 'LR'
            [junk,idx] = sort(real(lam),'descend');
        case 'SR'
            [junk,idx] = sort(real(lam));
        case 'LI'
            [junk,idx] = sort(imag(lam),'descend');
        case 'SI'
            [junk,idx] = sort(imag(lam));
        case 'LM'
            [junk,idx] = sort(abs(lam),'descend');
            % case 'SM' already converted to sigma = 0
        otherwise
            error('CHEBFUN:linop:eigs:sigma', 'Unidentified input ''sigma''.');
    end
end

% Delete infinite values. These can arise from rank deficiencies in the
% RHS matrix of the generalized eigenproblem.
idx( ~isfinite(lam(idx)) ) = [];

% Propose to keep these modes.
queue = 1:min(k,length(idx));
keeper = false(size(idx));
keeper(queue) = true;

%%
% Screen out spurious modes. These are dominated by high frequency for all
% values of N. (Known to arise for some formulations in generalized
% eigenproblems, specifically Orr-Sommerfeld.)

% Grab some indices
tenPercent = ceil(N/10); % We are the 10%
iif10 = 1:tenPercent;    % Indices of first 10%
ii90 = tenPercent:N;     % Indices of last 90%
ii10 = (N-tenPercent):N; % Indices of last 10%

% Check for high frequency energy (indicative of spurious eigenvalues) in
% each of the remaining valid eigenfunctions.
isFun = disc.source.isFunVariable;
while ~isempty(queue)
    j = queue(1);
    
    vc = mat2poly(disc,V(:,idx(j)));
    vc = vc(isFun);
    vcsq = 0;
    for i = 1:numel(vc)
        for q = 1:numel(vc{i})
            vcsq = vcsq + (vc{i}{q}.*conj(vc{i}{q}));
        end
    end
    vc = sqrt( flipud( sum(vcsq,2) ) ); 
      
    % Recipe: More than half of the energy in the last 90% of the Chebyshev
    % modes is in the highest 10% modes, and the energy of the last 90% is
    % not really small (1e-8) compared to the first 10% (i.e. is not noise).
    norm90 = norm(vc(ii90)); % Norm of last 90%
    norm10 = norm(vc(ii10)); % Norm of last 10%
    normfirst10 = norm(vc(iif10)); % Norm of first 10%
    if ( norm10 > 0.5*norm90 && norm90 > 1e-8*normfirst10 )
        keeper(j) = false;
        if queue(end) < length(idx)
            m = queue(end)+1;
            keeper(m) = true;  
            queue = [queue(:); m];
        end
    end
    queue(1) = [];
    
end

%%

% Return the keepers.
idx = idx( keeper );

end
