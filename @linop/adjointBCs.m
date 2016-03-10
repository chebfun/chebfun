function [Cstar, bcOpL, bcOpR, bcOpM] = adjointBCs(L, bcType)
%ADJOINTBCS   Compute the adjoint boundary conditions of a differential LINOP.
%   [Cstar, BCOPL, BCOPR, BCOPM] = ADJOINTBCS(L,BCTYPE), where L is a LINOP,
%   returns a set of functional blocks and function handles for the adjoint
%   boundary conditions of a LINOP L under the assumption that L only has
%   endpoint or periodic functional constraints.
%
% See also ?.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Get the domain and the value of the highest derivative:
dom = L.domain;
DOR = L.diffOrder;
dor = max(DOR);
[nout, nin] = size(L);

% initialize output
Cstar = L.constraint;
Bstar = [];

% initialize funs
funs = Cstar.functional;
nbcs = size(funs,1);
nadjbcs = 2*sum(dor)-nbcs;

% periodic case
if ( strcmp(bcType,'periodic') )
    bcOpL = [];
    bcOpR = [];
    bcOpM = 'periodic';

% bvp case
else

  % compute B matrix using constraints of L
  % construct chebmatrix of test functions
  fU = [];
  endValsL = [];
  endValsR = [];
  E = eye(nin);
  for ii = 1:nin
      % apply functionals to test basis
      bCol = [];
      z = zeros(1,2*dor(ii));
      U = chebpoly(0:2*dor(ii)-1,dom);
      for jj = 1:nin
           bCol = [ bCol; E(jj,ii)*U ];
      end
      fU = [ fU funs*bCol ];

      % evaluate test functions at endpoints
      Ul = U(dom(1),:); Ur = U(dom(2),:);
      for jj = 1:dor(ii)-1
          U = diff(U);
          Ul = [Ul;U(dom(1),:)]; Ur = [Ur;U(dom(2),:)];
      end
      endValsL = blkdiag( endValsL, Ul );
      endValsR = blkdiag( endValsR, Ur );
  end
  fU
  endVals = [ endValsL; endValsR ]

  % solve system to get B
  B = (endVals'\fU')'
  if ( rank(B) ~= nbcs )
    error('CHEBFUN:LINOP:adjoint:boundaryconditions', ...
    'Boundary condtions of L are not linearly independent.');
  end

  % attempt to simplify rows of B
  B = rref(B);
  B( abs(B-1) < 10*eps ) = 1;
  B( abs(B) < 10*eps )   = 0;
  B

  % label each boundary condition as 
  % left (0), right (1) or mixed (2)
  bcTypes = 2*ones(nbcs,1);
  for ii = 1:nbcs
    if ( max(abs(B(ii,1:sum(dor)))) == 0 )
      bcTypes(ii,:) = 1;
    elseif ( max(abs(B(ii,sum(dor)+1:end))) == 0 )
      bcTypes(ii,:) = 0;
    end
  end
  bcTypes

  % need the null space of B which we get via QR
  if ( sum(bcTypes == 2) )
    [qB,~] = qr(B');
    nulB = qB(:,nbcs+1:end);
  % if no mixed boundary conditions we ensure that the adjoint bcs 
  % are strictly left and right
  else
    nlbcs = sum(bcTypes == 0);
    [qBl,~] = qr(B(1:nlbcs,1:sum(dor))');
    nrbcs = nbcs - nlbcs;
    [qBr,~] = qr(B(nlbcs+1:end,sum(dor)+1:2*sum(dor))');
    nulB = zeros(2*sum(dor),nadjbcs);
    nulB(1:sum(dor),1:sum(dor)-nlbcs) = qBl(:,nlbcs+1:end);
    nulB(sum(dor)+1:2*sum(dor),sum(dor)-nlbcs+1:nadjbcs) = qBr(:,nrbcs+1:end);
  end
  nulB

  % compute complimentarity matrix
  DOR = L.diffOrder;
  compML = cell(nout,nin);
  compMR = cell(nout,nin);
  for ii = 1:nin
      for jj = 1:nout
          bloc = L.blocks{jj,ii};
          compML{jj,ii} = compmat( dom(1), toCoeff(bloc),...
                                  max(DOR(:,jj)), max(DOR(:,ii)) );
          compMR{jj,ii} = compmat( dom(2), toCoeff(bloc),...
                                  max(DOR(:,jj)), max(DOR(:,ii)) );
      end
  end
  compM = blkdiag( -cell2mat(compML), cell2mat(compMR) )

  % compute Bstar
  Bstar = nulB'*compM;
  Bstar	

  % attempt to simplify rows of Bstar
  Bstar = rref(Bstar);
  Bstar( abs(Bstar-1) < 10*eps ) = 1;
  Bstar( abs(Bstar) < 10*eps )   = 0;
  Bstar

  % label each boundary condition as 
  % left (0), right (1) or mixed (2)
  starTypes = 2*ones(nadjbcs,1);
  for ii = 1:nadjbcs
    if ( max(abs(Bstar(ii,1:sum(dor)))) == 0 )
      starTypes(ii,:) = 1;
    elseif ( max(abs(Bstar(ii,sum(dor)+1:end))) == 0 )
      starTypes(ii,:) = 0;
    end
  end
  starTypes

  % sort in ascending order
  [starTypes,ind] = sort(starTypes,'ascend');
  Bstar = Bstar(ind,:);
  Bstar

  % call bchandles
  [ bcOpL, bcOpR, bcOpM ] = bcHandles( Bstar, starTypes, dor );

  % compute Cstar by first creating a chebop
  Z = '[';
  for ii = 1:nin
      for jj = 1:nout
          Z = [Z,'+diff(v',int2str(jj),',',int2str(DOR(ii,jj)),')'];
      end
      Z = [Z,';'];
  end
  Z = [Z(1:end-1),']'];
  op = '@(x';
  for ii = 1:nout
      op = [op,',v',int2str(ii)];
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




function [bcOpL, bcOpR, bcOpM] = bcHandles(Bstar,starTypes,diffOrders)

% initialize ops
bcOpL = [];
bcOpR = [];
bcOpM = [];

% create cell array of dual variable names
nvars = length(diffOrders);
varNames = cell(1,nvars);
for ii = 1:nvars
    varNames{1,ii} = ['v',int2str(ii)];
end

% create cell array of diffnames
diffNames = cell(1,nvars);
for ii = 1:nvars
    set = {varNames{ii}};
    for jj = 1:diffOrders(ii)
        if ( jj == 1 )
            set{1,jj+1} = ['diff(',set{1,1},')'];
        else
            set{1,jj+1} = ['diff(',set{1,1},',',int2str(jj),')'];
        end
    end
    diffNames{ii} = set;
end

% convert each row of Bstar into a string regardless of starTypes
nbcs = size(Bstar,1);
nin = size(Bstar,2)/2;
Bstg = cell(nbcs,1);
for ii = 1:nbcs
    % left bcs
    for jj = 1:nvars
        stride = (jj-1)*diffOrders(max(jj-1,1));
        for kk = 1:diffOrders(jj)
            % add nonzero terms in row
            b = Bstar(ii,stride+kk);
            if ( b == 1 )
                % add plus sign
                if ( ~isempty(Bstg{ii,1}) )
                    Bstg{ii,1} = [Bstg{ii,1},'+'];
                end
                Bstg{ii,1} = [Bstg{ii,1},diffNames{jj}{kk}];
            elseif ( b ~= 0 )
                % add plus sign
                if ( ~isempty(Bstg{ii,1}) )
                    Bstg{ii,1} = [Bstg{ii,1},'+'];
                end
                bs = ['bl',int2str(ii),int2str(stride+kk)];
                eval([bs,'= b;']);
                Bstg{ii,1} = [Bstg{ii,1},bs,'*',diffNames{jj}{kk}];
            end
        end
    end
    % right bcs
    for jj = 1:nvars
        stride = nin+(jj-1)*diffOrders(max(jj-1,1));
        for kk = 1:diffOrders(jj)
            % add nonzero terms in row
            b = Bstar(ii,stride+kk);
            if ( b == 1 )
                % add plus sign
                if ( ~isempty(Bstg{ii,1}) )
                    Bstg{ii,1} = [Bstg{ii,1},'+'];
                end
                Bstg{ii,1} = [Bstg{ii,1},diffNames{jj}{kk}];
            elseif ( b ~= 0 )
                % add plus sign
                if ( ~isempty(Bstg{ii,1}) )
                    Bstg{ii,1} = [Bstg{ii,1},'+'];
                end
                bs = ['br',int2str(ii),int2str(stride+kk)];
                eval([bs,'= b;']);
                Bstg{ii,1} = [Bstg{ii,1},bs,'*',diffNames{jj}{kk}];
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
    bcOpL = ['@(',varNames{1,1}];
    for ii = 2:nvars
        bcOpL = [bcOpL,',',varNames{1,ii}];
    end
    bcOpL = [bcOpL,')'];
    % opening bracket if vectorized
    if ( length(indsL) > 1 )
        bcOpL = [bcOpL,'['];
    end
    % loop through indsL
    for ii = 1:length(indsL)
        bcOpL = [bcOpL,' ',Bstg{indsL(ii),1},';'];
    end 
    % closing bracket if vectorized
    if ( length(indsL) > 1 )
        bcOpL = [bcOpL(1:end-1),'];'];
    end
    % create function handle
    eval(['bcOpL = ',bcOpL]);
end
bcOpL

% bcOpR
if ( ~isempty(indsR) )
    % create variable list
    bcOpR = ['@(',varNames{1,1}];
    for ii = 2:nvars
        bcOpR = [bcOpR,',',varNames{1,ii}];
    end
    bcOpR = [bcOpR,')'];
    % opening bracket if vectorized
    if ( length(indsR) > 1 )
        bcOpR = [bcOpR,' ['];
    end
    % loop through indsR
    for ii = 1:length(indsR)
        bcOpR = [bcOpR,' ',Bstg{indsR(ii),1},';'];
    end 
    % closing bracket if vectorized
    if ( length(indsR) > 1 )
        bcOpR = [bcOpR(1:end-1),'];'];
    end
    % create function handle
    eval(['bcOpR = ',bcOpR]);
end
bcOpR

% bcOpM
if ( ~isempty(indsM) )
    % create variable list
    bcOpM = ['@(x,',varNames{1,1}];
    for ii = 2:nvars
        bcOpM = [bcOpM,',',varNames{1,ii}];
    end
    bcOpM = [bcOpM,')'];
    % opening bracket if vectorized
    if ( length(indsM) > 1 )
        bcOpM = [bcOpM,' ['];
    end
    % loop through indsM
    for ii = 1:length(indsM)
        bcOpM = [bcOpM,' ',Bstg{indsM(ii),1},';'];
    end 
    % closing bracket if vectorized
    if ( length(indsM) > 1 )
        bcOpM = [bcOpM(1:end-1),'];'];
    end
    % create function handle
    eval(['bcOpM = ',bcOpM]);
end
bcOpM

end




function A = compmat(x,coeffs,nrows,ncols)
%COMPMAT - routine for computing complementarity matrix.

  n = length(coeffs)-1;

  A = zeros(nrows,ncols);
  for ii=0:min(n-1,nrows-1)
    for jj=0:min(n-1-ii,ncols-1)
      for kk = 0:n
        if ( ii+jj <= kk-1 )
          A(ii+1,jj+1) = A(ii+1,jj+1) + (-1)^(kk-ii-1)*nchoosek(kk-ii-1,jj)*...
                         feval(diff(coeffs{n-kk+1},kk-ii-jj-1),x);
        end
      end
    end
  end

end
