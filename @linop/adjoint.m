function [Lstar, op, bcOpL, bcOpR, bcOpM] = adjoint(L, bcType)
%ADJOINT   Compute the adjoint of a LINOP.
%   ADJOINT(L), where L is a LINOP, returns the adjoint LINOP of L under
%   the assumption that L only has endpoint or periodic functional constraints.
%
%   ADJOINT(L, BCTYPE) allows for more general boundary conditions.
%
% See also ?.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Trivial case:
if ( L.diffOrder == 0 )
    Lstar = L;
    return
end

%% 
% Check to so if we can do compute an adjoint for this operator!
parseInputs(L, bcType);

%%
% Formal adjoint:
pref = chebfunpref();
if ( strcmpi(bcType, 'periodic') )
    pref.tech = @trigtech;
end 
[Lstar, op] = formalAdjoint(L, pref);

%%
% Create adjoint linop
if ( strcmp(bcType, 'periodic') )
    % Periodic bcs
    Lstar.constraint = L.constraint;
    bcOpL = [];
    bcOpR = [];
    bcOpM = 'periodic';
else
    % Adjoint boundary conditions
    [constraint, bcOpL, bcOpR, bcOpM] = adjointBCs(L, bcType);
    Lstar.constraint = constraint;
end

end

function parseInputs(L, bcType)

if ( L.diffOrder < 0 )
    % Check for integral operators:
    error('CHEBFUN:LINOP:adjoint:difforder', ...
    'ADJOINT doesn''t support integral operators for the moment.')
elseif ( size(L.domain, 2) > 2 )
    % [TODO]: Support piecewise domains.
    error('CHEBFUN:LINOP:adjoint:domain', ...
        'ADJOINT doesn''t support piecewise domains for the moment.');
elseif ( ~any(strcmp(bcType, {'periodic', 'bvp'})) )
    error('CHEBFUN:LINOP:adjoint:boundaryconditions', ...
    'ADJOINT doesn''t support this type of boundary conditions for the moment.');
elseif ( all(size(L.blocks)) ~= 1 )
    % [TODO]: Support systems of equations.
    error('CHEBFUN:LINOP:adjoint:notscalar', ...
    'ADJOINT only support square scalar systems..');
end

end

function [L, op] = formalAdjoint(L, pref)
%FORMALADJOINT   Computes the formal adjoint of a LINOP.
%   L = FORMALADJOINT(L, BCTYPE) returns the formal adjoint of a LINOP L.
%
% See also adjointBCs.

% Get the domain and the value of the highest derivative:
dom = L.domain;
n = L.diffOrder;

% Get the coefficients:
coeffs = toCoeff(L.blocks{1}, pref);

% Compute the coefficients of the adjoint:
adjCoeffs = 0*coeffs;
for k = 0:n
    for l = 0:k
        adjCoeffs(n+1-l) = adjCoeffs{n+1-l} + ...
            (-1)^k*nchoosek(k,l)*conj(diff(coeffs{n+1-k}, k-l));
    end
end

% Construct a LINOP from these new coefficients:
L = 0;
M = @(f) operatorBlock.mult(f, dom);
D = @(k) operatorBlock.diff(dom, k);
a = [];
for k = 0:n
    varname = genvarname(['a', int2str(k)]);
    eval([varname, '= adjCoeffs{n+1-k};']);
    L = L + M(varname) * D(k);
end
L = linop(L);

% set op
op = 'a0*u';
for k = 1:n
  op = ['a',int2str(k),'*diff(u,',int2str(k),') + ',op];
end
eval(['op = @(x,u) ',vectorize(op),';']);

end





function [constraint, bcOpL, bcOpR, bcOpM] = adjointBCs(L, bcType)
%ADJOINTBCS   Computes the adjoint boundary conditions of a LINOP.
%
% See also FORMALADJOINT.

% initialize output
constraint = L.constraint;
Bstar = [];

% return if bcType is periodic
if ( strcmp(bcType,'periodic') )
    return
end

% if bcType ~= 'bvp' throw error
if ( ~strcmp(bcType, 'bvp') )
    error('CHEBFUN:LINOP:adjoint:boundaryconditions', ...
    'ADJOINT doesn''t support this type of boundary conditions for the moment.');
end

% Get the domain and the value of the highest derivative:
dom = L.domain;
n = L.diffOrder;

% Get the coefficients:
coeffs = toCoeff(L.blocks{1});

% compute complementary matrix
A11 = compmat(dom(1),coeffs);
A22 = compmat(dom(2),coeffs);
M = zeros(2*n);
M(1:n,1:n) = -A11;
M(n+1:end,n+1:end) = A22;

% compute B matrix using constraints of L
% initialize funs
funs = constraint.functional;
nbcs = size(funs,1);
nadjbcs = 2*n-nbcs;

% vector of chebpolys and derivatives
V = chebpoly(0:2*n-1,dom); F = funs*V;
vl = V(dom(1),:); vr = V(dom(2),:);
for kk = 1:n-1
  V = diff(V);
  vl = [vl;V(dom(1),:)]; vr = [vr;V(dom(2),:)];
end
v = [vl;vr];

% solve system to get B
B = (v'\F')';
if ( rank(B) ~= nbcs )
    error('CHEBFUN:LINOP:adjoint:boundaryconditions', ...
    'Boundary condtions of L are not linearly independent.');
end
[~,ind] = max(abs(B),[],2);
rmax = diag(B(1:end,ind));
B = diag(1./rmax)*B;
B( abs(B) < 10*eps ) = 0;

% label each boundary condition as 
% left (0), right (1) or mixed (2)
bcTypes = 2*ones(nbcs,1);
for ii = 1:nbcs
  if ( max(abs(B(ii,1:n))) == 0 )
    bcTypes(ii,:) = 1;
  elseif ( max(abs(B(ii,n+1:end))) == 0 )
    bcTypes(ii,:) = 0;
  end
end

% need the null space of B which we get via QR
if ( sum(bcTypes == 2) )
  [qB,~] = qr(B');
  nulB = qB(:,nbcs+1:end);
% if no mixed boundary conditions we ensure that the adjoint bcs 
% are strictly left and right
else
  nlbcs = sum(bcTypes == 0);
  [qBl,~] = qr(B(1:nlbcs,1:n)');
  nrbcs = nbcs - nlbcs;
  [qBr,~] = qr(B(nlbcs+1:end,n+1:2*n)');
  nulB = zeros(2*n,nadjbcs);
  nulB(1:n,1:n-nlbcs) = qBl(:,nlbcs+1:end);
  nulB(n+1:2*n,n-nlbcs+1:nadjbcs) = qBr(:,nrbcs+1:end);
end

% compute Bstar
Bstar = nulB'*M;
[~,ind] = max(abs(Bstar),[],2);
rmax = diag(Bstar(1:end,ind));
Bstar = diag(1./rmax)*Bstar;
Bstar( abs(Bstar) < 10*eps ) = 0;

% label each boundary condition as 
% left (0), right (1) or mixed (2)
starTypes = 2*ones(nadjbcs,1);
for ii = 1:nadjbcs
  if ( max(abs(Bstar(ii,1:n))) == 0 )
    starTypes(ii,:) = 1;
  elseif ( max(abs(Bstar(ii,n+1:end))) == 0 )
    starTypes(ii,:) = 0;
  end
end

% sort in ascending order
[starTypes,ind] = sort(starTypes,'ascend');
Bstar = Bstar(ind,:);

% create star functionals
FL = functionalBlock.feval(dom(1), dom);
FR = functionalBlock.feval(dom(2), dom);
D = @(k) operatorBlock.diff(dom, k);
starFuns = cell(nadjbcs,1);
for ii = 1:nadjbcs
  starFuns{ii} = 0*FL;
  % left bcs
  if ( starTypes(ii) == 0 || starTypes(ii) == 2 )
    for jj = 1:n
      starFuns{ii} = starFuns{ii} + Bstar(ii,jj)*FL*D(jj-1);
    end
  end
  % right bcs
  if ( starTypes(ii) == 1 || starTypes(ii) == 2 )
    for jj = 1:n
      starFuns{ii} = starFuns{ii} + Bstar(ii,n+jj)*FR*D(jj-1);
    end
  end
end

% create output constraints
constraint.functional = chebmatrix(starFuns);
constraint.values = zeros(nadjbcs,1);

% create function handles for left, right and mixed bcs
inds = (1:nadjbcs)';
indsL = inds( starTypes == 0 )';
indsR = inds( starTypes == 1 )';
indsM = inds( starTypes == 2 )';
bcOpL = [];
bcOpR = [];
bcOpM = [];

% left
if ( ~isempty(indsL) )
bcOpL = '@(u) [';
for ii = indsL
  for jj = 1:n
    if ( Bstar(ii,jj) ~= 0 )
      bcOpL = [bcOpL,' + ',num2str(Bstar(ii,jj),'%1.15e'),...
               '*diff(u,',int2str(jj-1),')'];
    end
  end
  bcOpL = [bcOpL,'; '];
end
bcOpL = [bcOpL,']'];
bcOpL = str2func(vectorize(bcOpL));
end

% right
if ( ~isempty(indsR) )
bcOpR = '@(u) [';
for ii = indsR
  for jj = 1:n
    if ( Bstar(ii,n+jj) ~= 0 )
      bcOpR = [bcOpR,' + ',num2str(Bstar(ii,n+jj),'%1.15e'),...
               '*diff(u,',int2str(jj-1),')'];
    end
  end
  bcOpR = [bcOpR,'; '];
end
bcOpR = [bcOpR,']'];
bcOpR = str2func(vectorize(bcOpR));
end

% mixed
if ( ~isempty(indsM) )
bcOpM = '@(x,u) [';
for ii = indsM
  for jj = 1:n
    if ( Bstar(ii,jj) ~= 0 )
      bcOpM = [bcOpM,' + ',num2str(Bstar(ii,jj),'%1.15e'),...
               '*deriv(u,',num2str(dom(1),'%1.15e'),...
               ',',int2str(jj-1),')'];
    end
  end
  for jj = 1:n
    if ( Bstar(ii,n+jj) ~= 0 )
      bcOpM = [bcOpM,' + ',num2str(Bstar(ii,n+jj),'%1.15e'),...
               '*deriv(u,',num2str(dom(2),'%1.15e'),...
               ',',int2str(jj-1),')'];
    end
  end
  bcOpM = [bcOpM,'; '];
end
bcOpM = [bcOpM,']'];
bcOpM = str2func(vectorize(bcOpM));
end

end




function A = compmat(x,coeffs)
%COMPMAT - routine for computing complementarity matrix.

  n = length(coeffs)-1;

  A = zeros(n);
  for ii=0:n-1
    for jj=0:n-1-ii
      for kk = 0:n
        if ( ii+jj <= kk-1 )
          A(ii+1,jj+1) = A(ii+1,jj+1) + (-1)^(kk-ii-1)*nchoosek(kk-ii-1,jj)*...
                         feval(diff(coeffs{n-kk+1},kk-ii-jj-1),x);
        end
      end
    end
  end

  nrmA = norm(A);
  if ( nrmA ~= 0 )
    A = A/nrmA;
  end

end
