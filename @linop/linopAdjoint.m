function [Lstar, op, bcOpL, bcOpR, bcOpM] = linopAdjoint(L, bcType)
%LINOPADJOINT   Compute the adjoint of a LINOP.
%   [LSTAR, OP, BCOPL, BCOPR, BCOPM] = LINOPADJOINT(L, BCTYPE) computes the 
%   adjoint of the LINOP L. LINOPADJOINT requires that L represents a linear 
%   differential operator with either endpoint or periodic boundary conditions.
%   If L represents a system of differential equations then the highest order
%   derivative in each variable must be the same. Integral operators and exotic
%   functional constraints are not supported.
%
%   The output is a LINOP LSTAR that represents the adjoint differential
%   operator and adjoint constraints, as well as four function handles that can
%   be used to construct a CHEBOP. OP is a function handle for the differential
%   operator and BCOPL, BCOPR and BCOPM are function handles for the left,
%   right and mixed boundary conditions respectively. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check for 2 inputs:
if ( nargin < 2 )
    error('CHEBFUN:LINOP:linopAdjoint:inputs', ...
       'LINOPADJOINT requires two inputs.')
end

% Check to see if we can compute an adjoint for this operator!
parseInputs(L, bcType);

% Formal adjoint:
pref = chebfunpref();
if ( strcmpi(bcType, 'periodic') )
    pref.tech = @trigtech;
end 
[Lstar, op] = adjointFormal(L, pref);

% Trivial case:
if ( max(max(L.diffOrder)) == 0 )
    bcOpL = [];
    bcOpR = [];
    bcOpM = [];
    return
end

% Adjoint boundary conditions:
[Cstar, bcOpL, bcOpR, bcOpM] = adjointBCs(L, bcType);
Lstar.constraint = Cstar;

end

function parseInputs(L, bcType)
%PARSEINPUTS   Function to parse inputs and catch errors.

[m,n] = size(L);

if ( min(min(L.diffOrder)) < 0 )
    % Check for integral operators:
    error('CHEBFUN:LINOP:adjoint:difforder', ...
        'LINOPADJOINT doesn''t support integral operators for the moment.')
elseif ( size(L.domain, 2) > 2 )
    % [TODO]: Support piecewise domains.
    error('CHEBFUN:LINOP:adjoint:domain', ...
        'LINOPADJOINT doesn''t support piecewise domains for the moment.');
elseif ( ~any(strcmp(bcType, {'periodic', 'bvp', 'ivp', 'fvp'})) )
    error('CHEBFUN:LINOP:adjoint:boundaryconditions', ...
        ['LINOPADJOINT doesn''t support this type of boundary conditions', ...
         'for the moment.']);
elseif ( m ~= n )
    % [TODO]: Support nonsquare systems.
    error('CHEBFUN:LINOP:adjoint:systems', ...
        'LINOPADJOINT doesn''t support nonsquare systems at the moment.');
elseif ( n > 1 && any(diff(max(L.diffOrder)) ~= 0) )
    % [TODO]: Support block systems with arbitrary differential order in 
    %         each variable.
    error('CHEBFUN:LINOP:adjoint:systems', ...
        ['LINOPADJOINT doesn''t support systems with arbitrary ', ...
         'differential order in each variable at the moment.']);
end

end

function [Lstar, op] = adjointFormal(L, pref)
%ADJOINTFORMAL   Compute the formal adjoint of a LINOP.
%   [LSTAR, OP] = ADJOINTFORMAL(L, PREF), where L is a LINOP, returns the
%   formal adjoint LINOP of L, LSTAR, under the condition that L is a linear
%   differential operator, i.e. 
%
%       L*u = a_k*diff(u,k) + ... + a_0*u.
%
%   If L represents a system of differential equations then the differential
%   order in each variable must be the same. Integral operators are not
%   supported. 
%
%   The outputs are a LINOP LSTAR that represents the formal adjoint and a
%   function handle OP that can be used to construct a CHEBOP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check for 2 inputs:
if ( nargin < 2 )
    error('CHEBFUN:LINOP:adjointFormal:inputs', ...
        'LINOPADJOINT requires two inputs.')
end

% Initialize Lstar and op to have correct dimensions:
[nrows, ncols] = size(L);
op = mat2cell(zeros(ncols, nrows), ones(1, ncols), ones(1, nrows));
Lstar = linop(op);

% If L.blocks{ii,jj} is a chebfun convert to operator block:
for ii = 1:nrows
    for jj = 1:ncols
        if ( isa(L.blocks{ii,jj}, 'chebfun') )
            b = L.blocks{ii,jj};
            L.blocks{ii,jj} = operatorBlock.mult(b, b.domain);
        end
    end
end

% Get the domain and the value of the highest derivative:
dom = L.domain;
Lstar.domain = dom;

% Loop through blocks:
op = ';';
argstr = '@(x';
for ii = 1:ncols
    for jj = 1:nrows
        
        % Get local block:
        B = L.blocks{jj,ii};

        % Get the coefficients:
        coeffs = toCoeff(B, pref);
        dor = length(coeffs)-1;

        % Compute the coefficients of the adjoint:
        adjCoeffs = 0*coeffs;
        for k = 0:dor
            for l = 0:k
                adjCoeffs(dor+1-l) = adjCoeffs{dor+1-l} + ...
                    (-1)^k*nchoosek(k,l)*conj(diff(coeffs{dor+1-k}, k-l));
            end
        end

        % Update the blocks of Lstar using adjoint coefficients:
        Lstar.blocks{ii,jj} = 0;
        M = @(f) operatorBlock.mult(f, dom);
        D = @(k) operatorBlock.diff(dom, k);
        for k = 0:dor
            varname = ['a', int2str(ii), int2str(jj), '_', int2str(k)];
            eval([varname, '= adjCoeffs{dor+1-k};']);
            Lstar.blocks{ii,jj} = Lstar.blocks{ii,jj} + ...
				  M(adjCoeffs{dor+1-k}) * D(k);
        end

        % Update argstr:
        if nrows == 1
            vstr = 'v';
        else
            vstr = ['v', int2str(jj)];
        end
        if ( ii == 1 )
            argstr = [argstr, ',', vstr];
        end

        % Set Bop:
        Bop = [];
        astr = ['a', int2str(ii), int2str(jj)];
        for k = 0:dor
            
            acur = [astr, '_', int2str(k)];
            
            % Check acur for 0, 1, or -1:
            acheb = eval(acur);
            if ( norm(acheb) == 0 && dor == 0 )
                acur = '+0*';
            elseif ( norm(acheb) == 0 && dor > 0 )
                acur = [];
            elseif ( length(acheb) == 1 && acheb(0) == 1 )
                acur = '+';
            elseif ( length(acheb) == 1 && acheb(0) == -1 )
                acur = '-';
            else 
                acur = ['+', acur, '.*'];
            end
            
            % Use acur to update op string:
            if ( ~isempty(acur) )
                if ( k == 0 )
                    Bop = [acur, vstr, Bop];
                elseif ( k == 1 )
                    Bop = [acur, 'diff(', vstr, ')', Bop];
                else
                    Bop = [acur, 'diff(', vstr, ',', int2str(k), ')', Bop];
                end
            end
        end
 
        % Update op:
        if ( ~isempty(Bop) && strcmp(op(end), ';') )
            if ( strcmp(Bop(1), '+') )
                Bop = [' ', Bop(2:end)];
            elseif ( strcmp(Bop(1), '-') )
                Bop = [' ', Bop];
            end
        end
        op = [op, Bop];

    end

    % Update op:
    if ( ii < ncols )
        op = [op, ';'];
    end

end 

% Set op:
if ( ncols > 1 )
    op = [argstr, ') [', op(2:end), ' ];'];
else
    op = [argstr, ') ', op(2:end), ';'];
end
eval(['op = ', op]);

end

function [Cstar, bcOpL, bcOpR, bcOpM] = adjointBCs(L, bcType)
%ADJOINTBCS   Computes the adjoint boundary conditions of a LINOP.
%   [CSTAR, BCOPL, BCOPR, BCOPM] = ADJOINTBCS(L, BCTYPE), where L is a LINOP,
%   returns a set of functional blocks and function handles for the adjoint
%   boundary conditions of a LINOP L under the assumption that L only has
%   endpoint or periodic functional constraints. If L represents a system
%   of differential equations then the differential order in each variable
%   must be the same. 
%
%   The output is a set of functional blocks CSTAR that can be used to
%   construct a LINOP and function handles BCOPL, BCOPR and BCOPM that
%   represent the left, right and mixed boundary conditions respectively and
%   can be used to construct a CHEBOP

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Get the domain and the value of the highest derivative:
dom = L.domain;
DOR = L.diffOrder;
dor = max(DOR);
[nout, nin] = size(L);

% initialize output
Cstar = L.constraint;
bcOpL = [];
bcOpR = [];
bcOpM = [];

% initialize funs
funs = Cstar.functional;
nbcs = size(funs,1);
nadjbcs = 2*sum(dor)-nbcs;

% periodic case
if ( strcmp(bcType,'periodic') )
    bcOpM = 'periodic';

% bvp case
else

  % compute B matrix using constraints of L
  % construct chebmatrix of test functions
  fU = [];
  endValsL = [];
  endValsR = [];
  E = eye(nin);
  warning('off','all') % is there a better way to supress warnings
  for ii = 1:nin
      % apply functionals to test basis
      bCol = [];
      z = zeros(1,2*dor(ii));
      U = chebpoly(0:2*dor(ii)-1,dom);
      for jj = 1:nin
           bCol = [ bCol; E(jj,ii)*U ];
      end
      if ( isempty(funs) )
          fU = [ fU 0*feval(bCol,dom(1)) ];
      else
          fU = [ fU funs*bCol ];
      end

      % evaluate test functions at endpoints
      Ul = U(dom(1),:); 
      Ur = U(dom(2),:);
      for jj = 1:dor(ii)-1
          U = diff(U);
          Ul = [Ul; U(dom(1),:)]; 
          Ur = [Ur; U(dom(2),:)];
      end
      endValsL = blkdiag(endValsL, Ul);
      endValsR = blkdiag(endValsR, Ur);
  end
  warning('on','all') % is there a better way to supress warnings
  endVals = [ endValsL; endValsR ];

  % solve system to get B
  B = (endVals'\fU')';
  if ( rank(B) ~= nbcs )
    error('CHEBFUN:LINOP:adjoint:boundaryconditions', ...
    'Boundary condtions of L are not linearly independent.');
  end

  % attempt to simplify rows of B
  B = rref(B);
  B( abs(B-1) < 10*eps ) = 1;
  B( abs(B) < 10*eps ) = 0;

  % label each boundary condition as left (0), right (1) or mixed (2)
  bcTypes = 2*ones(nbcs,1);
  for ii = 1:nbcs
    if ( max(abs(B(ii,1:sum(dor)))) == 0 )
      bcTypes(ii,:) = 1;
    elseif ( max(abs(B(ii,sum(dor)+1:end))) == 0 )
      bcTypes(ii,:) = 0;
    end
  end

  % need the null space of B which we get via QR
  [qB,~] = qr(B');
  nulB = qB(:,nbcs+1:end)';
  nulB = rref(nulB);

  % compute complimentarity matrix
  compML = cell(nout, nin);
  compMR = cell(nout, nin);
  for ii = 1:nout
      for jj = 1:nin
          bloc = L.blocks{ii,jj};
          compML{ii,jj} = compmat(dom(1), toCoeff(bloc));
          compMR{ii,jj} = compmat(dom(2), toCoeff(bloc));
      end
  end
  compM = blkdiag(-cell2mat(compML), cell2mat(compMR));

  % compute Bstar
  Bstar = nulB*compM;

  % attempt to simplify rows of Bstar
  Bstar = rref(Bstar);
  Bstar( abs(Bstar-1) < 10*eps ) = 1;
  Bstar( abs(Bstar) < 10*eps ) = 0;

  % label each boundary condition as 
  % left (0), right (1) or mixed (2)
  starTypes = 2*ones(nadjbcs,1);
  for ii = 1:nadjbcs
    if ( max(abs(Bstar(ii, 1:sum(dor)))) == 0 )
      starTypes(ii,:) = 1;
    elseif ( max(abs(Bstar(ii, sum(dor)+1:end))) == 0 )
      starTypes(ii,:) = 0;
    end
  end

  % sort in ascending order
  [starTypes, ind] = sort(starTypes, 'ascend');
  Bstar = Bstar(ind,:);

  % call bchandles
  [ bcOpL, bcOpR, bcOpM ] = bcHandles(Bstar, starTypes, dor);

  % compute Cstar by first creating a chebop
  order = max(dor, [], 2);
  Z = '[';
  for ii = 1:nin
      for jj = 1:nout
          Z = [Z, '+diff(v', int2str(jj), ',', int2str(order), ')'];
      end
      Z = [Z,';'];
  end
  Z = [Z(1:end-1), ']'];
  op = '@(x';
  for ii = 1:nout
      op = [op, ',v', int2str(ii)];
  end
  op = [op,') ',Z,';'];
  N = chebop(dom);
  eval(['N.op = ',op]);
  N.lbc = bcOpL;
  N.rbc = bcOpR;
  N.bc = bcOpM;
  NL = linearize(N);
  Cstar = NL.constraint;

end 

end

function [bcOpL, bcOpR, bcOpM] = bcHandles(Bstar, starTypes, diffOrders)
%BCHANDLES   Construct function handles using Bstar matrix.

% initialize ops
bcOpL = [];
bcOpR = [];
bcOpM = [];

% variables for left and right endpoint
eval('lpt = domain(1);');
eval('rpt = domain(2);');

% create cell array of dual variable names
nvars = length(diffOrders);
varNames = cell(1, nvars);
if nvars == 1
    varNames = {'v'};
else
    for ii = 1:nvars
        varNames{1,ii} = ['v', int2str(ii)];
    end
end

% create cell array of diffnames
diffNames = cell(1, nvars);
for ii = 1:nvars
    set = {varNames{ii}};
    for jj = 1:diffOrders(ii)
        if ( jj == 1 )
            set{1,jj+1} = ['diff(', set{1,1}, ')'];
        else
            set{1,jj+1} = ['diff(', set{1,1}, ',', int2str(jj), ')'];
        end
    end
    diffNames{ii} = set;
end

% create cell array of diffnames at left end point
diffNamesL = cell(1, nvars);
for ii = 1:nvars
    set = {['feval(',varNames{ii},',lpt)']};
    for jj = 1:diffOrders(ii)
        if ( jj == 1 )
            set{1,jj+1} = ['feval(diff(', varNames{ii}, '),lpt)'];
        else
            set{1,jj+1} = ['feval(diff(', varNames{ii}, ',', int2str(jj), ...
                '),lpt)'];
        end
    end
    diffNamesL{ii} = set;
end

% create cell array of diffnames at right end point
diffNamesR = cell(1, nvars);
for ii = 1:nvars
    set = {['feval(', varNames{ii}, ',rpt)']};
    for jj = 1:diffOrders(ii)
        if ( jj == 1 )
            set{1,jj+1} = ['feval(diff(', varNames{ii}, '),rpt)'];
        else
            set{1,jj+1} = ['feval(diff(', varNames{ii}, ',', int2str(jj), ...
                '),rpt)'];
        end
    end
    diffNamesR{ii} = set;
end

% convert each row of Bstar into a string regardless of starTypes
nbcs = size(Bstar, 1);
nin = size(Bstar, 2)/2;
Bstg = cell(nbcs, 1);
for ii = 1:nbcs
    % choose list of diffnames
    if ( starTypes(ii) == 2 )
       dnames = diffNamesL;
    else
       dnames = diffNames;
    end
    % left bcs
    for jj = 1:nvars
        stride = (jj-1)*diffOrders(max(jj-1,1));
        for kk = 1:diffOrders(jj)
            % add nonzero terms in row
            b = Bstar(ii,stride+kk);
            if ( b == 1 )
                % add plus sign
                if ( ~isempty(Bstg{ii,1}) )
                    Bstg{ii,1} = [Bstg{ii,1}, '+'];
                end
                Bstg{ii,1} = [Bstg{ii,1}, dnames{jj}{kk}];
            elseif ( b ~= 0 )
                % add plus sign
                if ( ~isempty(Bstg{ii,1}) )
                    Bstg{ii,1} = [Bstg{ii,1}, '+'];
                end
                bs = ['bl',int2str(ii),int2str(stride+kk)];
                eval([bs,'= b;']);
                Bstg{ii,1} = [Bstg{ii,1}, bs, '*', dnames{jj}{kk}];
            end
        end
    end
    % choose list of diffnames
    if ( starTypes(ii) == 2 )
       dnames = diffNamesR;
    else
       dnames = diffNames;
    end
    % right bcs
    for jj = 1:nvars
        stride = nin+(jj-1)*diffOrders(max(jj-1, 1));
        for kk = 1:diffOrders(jj)
            % add nonzero terms in row
            b = Bstar(ii,stride+kk);
            if ( b == 1 )
                % add plus sign
                if ( ~isempty(Bstg{ii,1}) )
                    Bstg{ii,1} = [Bstg{ii,1}, '+'];
                end
                Bstg{ii,1} = [Bstg{ii,1}, dnames{jj}{kk}];
            elseif ( b ~= 0 )
                % add plus sign
                if ( ~isempty(Bstg{ii,1}) )
                    Bstg{ii,1} = [Bstg{ii,1}, '+'];
                end
                bs = ['br', int2str(ii), int2str(stride+kk)];
                eval([bs,'= b;']);
                Bstg{ii,1} = [Bstg{ii,1}, bs, '*', dnames{jj}{kk}];
            end
        end
    end
end

% create function handles based on starTypes
inds = (1:nbcs)';
indsL = inds(starTypes == 0);
indsR = inds(starTypes == 1);
indsM = inds(starTypes == 2);

% bcOpL
if ( ~isempty(indsL) )
    % create variable list
    bcOpL = ['@(', varNames{1,1}];
    for ii = 2:nvars
        bcOpL = [bcOpL, ',', varNames{1,ii}];
    end
    bcOpL = [bcOpL, ')'];
    % opening bracket if vectorized
    if ( length(indsL) > 1 )
        bcOpL = [bcOpL, '['];
    end
    % loop through indsL
    for ii = 1:length(indsL)
        bcOpL = [bcOpL, ' ', Bstg{indsL(ii),1}, ';'];
    end 
    % closing bracket if vectorized
    if ( length(indsL) > 1 )
        bcOpL = [bcOpL(1:end-1), '];'];
    end
    % create function handle
    eval(['bcOpL = ', bcOpL]);
end

% bcOpR
if ( ~isempty(indsR) )
    % create variable list
    bcOpR = ['@(', varNames{1,1}];
    for ii = 2:nvars
        bcOpR = [bcOpR, ',', varNames{1,ii}];
    end
    bcOpR = [bcOpR,')'];
    % opening bracket if vectorized
    if ( length(indsR) > 1 )
        bcOpR = [bcOpR,' ['];
    end
    % loop through indsR
    for ii = 1:length(indsR)
        bcOpR = [bcOpR, ' ', Bstg{indsR(ii),1}, ';'];
    end 
    % closing bracket if vectorized
    if ( length(indsR) > 1 )
        bcOpR = [bcOpR(1:end-1), '];'];
    end
    % create function handle
    eval(['bcOpR = ', bcOpR]);
end

% bcOpM
if ( ~isempty(indsM) )
    % create variable list
    bcOpM = ['@(x,', varNames{1,1}];
    for ii = 2:nvars
        bcOpM = [bcOpM, ',', varNames{1,ii}];
    end
    bcOpM = [bcOpM, ')'];
    % opening bracket if vectorized
    if ( length(indsM) > 1 )
        bcOpM = [bcOpM, ' ['];
    end
    % loop through indsM
    for ii = 1:length(indsM)
        bcOpM = [bcOpM, ' ', Bstg{indsM(ii),1}, ';'];
    end 
    % closing bracket if vectorized
    if ( length(indsM) > 1 )
        bcOpM = [bcOpM(1:end-1), '];'];
    end
    % create function handle
    eval(['bcOpM = ', bcOpM]);
end

end

function A = compmat(x,coeffs)
%COMPMAT   Routine for computing complementarity matrix.

  n = length(coeffs)-1;
  if n == 0
      A = 0;
      return
  end

  A = zeros(n);
  for kk = n:-1:1
  % loop over diff orders
    cf = coeffs{n-kk+1};
    for ii = 0:kk-1
    % loop through rows
      for jj = 0:kk-1-ii
      % loop through columns
        if ( ii+jj <= kk-1 )
	  A(ii+1,jj+1) = A(ii+1,jj+1) + (-1)^(kk-ii-1)*nchoosek(kk-ii-1,jj)*...
                         feval(diff(cf,kk-ii-jj-1),x);
        end
      end
    end
  end

end
