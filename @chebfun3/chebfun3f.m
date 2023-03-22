function f = chebfun3f(f, op, pref, dom, vectorize)
%CHEBFUN3F  Default CHEBFUN3 constructor.
%   Given a function OP of three variables, this code represents it as a
%   CHEBFUN3 object. A CHEBFUN3 object is a low-rank representation
%   expressing a function as a trilinear product of a discrete core tensor
%   and three quasimatrices consisting of univariate functions.
%
%   The CHEBFUN3F algorithm for constructing a CHEBFUN3 object has the
%   potential to require fewer function evaluations compared to the default
%   constructor.
%
%   The CHEBFUN3F constructor has three phases:
%
%   PHASE 1: The first phase attempts to identify fibers on a coarse grid
%   which approximate the span of the tensor of function evaluations in the
%   corresponding mode.
%
%   PHASE 2: The second phase attempts to resolve the corresponding
%   fibers by sampling along the fibers until the Chebyshev coefficients of
%   the fibers fall below machine epsilon.
%
%   PHASE 3: The third phase attempts to construct factor matrices and a
%   core tensor approximating the evaluation tensor on a fine grid. The
%   columns of the factor matrices are then converted into CHEBFUN objects.
%
%  The Chebfun3f algorithm is fully described in:
%   [1] S. Dolgov, D. Kressner, C. Stroessner, Functional Tucker
%   approximation using Chebyshev interpolation, SIAM J. Sci. Comput., 43
%   (2021), A2190-A2210.
%
% See also CHEBFUN3, CHEBFUN2, CHEBFUN3T and CHEBFUN3V.

% Copyright 2023 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set preferences:
tech             = pref.tech();
prefStruct       = pref.cheb3Prefs;
tpref            = tech.techPref;
grid             = tpref.minSamples;
maxSample        = tpref.maxLength;     % max polynomialDeg (not implemented)
maxSamplePhase1  = 363;                 % max coarseResolution (not implemented)
maxRank          = prefStruct.maxRank;  % max rank (not implemented)
maxRestarts      = 10;                  % max restarts after the sample test
pseudoLevel      = prefStruct.chebfun3eps;
passSampleTest   = prefStruct.sampleTest;

% Initialize:
n                = [grid, grid, grid];  % coarseResolution
m                = n;                   % fineResolution
r                = [6,6,6];             % rank
absTol              = pseudoLevel;
chebX            = @(i,n) dom(1) + ((-cos((i-1).*pi/(n-1))) + 1)*(dom(2)-dom(1))/2 ;
chebY            = @(i,n) dom(3) + ((-cos((i-1).*pi/(n-1))) + 1)*(dom(4)-dom(3))/2 ;
chebZ            = @(i,n) dom(5) + ((-cos((i-1).*pi/(n-1))) + 1)*(dom(6)-dom(5))/2 ;
reffun           = @(n) floor(sqrt(2)^(floor(2*log2(n)) + 1)) + 1;
restarts         = 0;
f.domain         = dom;

if isa(tech,'trigtech')
    %use equispace points instead
    chebX = @(i,n) dom(1) +  (i-1)/(n) *(dom(2)-dom(1));
    chebY = @(i,n) dom(3) +  (i-1)/(n) *(dom(4)-dom(3));
    chebZ = @(i,n) dom(5) +  (i-1)/(n) *(dom(6)-dom(5));
end

%% Main Loop
happy = 0;
while ~happy
    
    %% Phase 1
    happyPhase1 = 0;
    while ( ~happyPhase1 )
        J = initializeIndexRandomly(r(2), n(2), restarts);
        K = initializeIndexRandomly(r(3), n(3), restarts);
        
        % Handle to evaluate tensor entries of T_c
        ff = @(i,j,k) op(chebX(i,n(1)),chebY(j,n(2)),chebZ(k,n(3)));
        
        for iterations = 1:2
            happyPhase1 = 1;
            
            % ACA 1
            T1 = evalTensor(1:n(1),J,K,ff,vectorize);
            
            % Does the function blow up or evaluate to NaN?:
            if ( isinf(max(abs(T1(:)))) )
                error('CHEBFUN:CHEBFUN3:chebfun3f:inf', ...
                    'Function returned INF when evaluated');
            elseif ( any(isnan(T1(:))) )
                error('CHEBFUN:CHEBFUN3:chebfun3f:nan', ...
                    'Function returned NaN when evaluated');
            end
  
            T1 = reshape(T1,n(1),r(2)*r(3));
            [~, absTol] = getTol(T1, pseudoLevel, absTol,dom(2)-dom(1),tech);
            [Uc, ~, ~, I,I2] = ACA(T1, absTol, n(1));
            r(1) = size(I,2);
            JT1 = J;
            KT1 = K;
            
            % ACA 2
            T2 = evalTensor(I,1:n(2),K,ff,vectorize);
            T2 = reshape(permute(T2,[2,1,3]),n(2),r(1)*r(3));
            [~, absTol] = getTol(T2, pseudoLevel, absTol,dom(4)-dom(3),tech);
            [Vc, ~, ~, J, J2] = ACA(T2, absTol, n(2));
            r(2) = size(J,2);
            KT2 = K;
            IT2 = I;
            
            % ACA 3
            T3 = evalTensor(I,J,1:n(3), ff,vectorize);
            T3 = reshape(permute(T3,[3,1,2]),n(3),r(1)*r(2));
            [relTol, absTol] = getTol(T3, pseudoLevel, absTol,dom(6)-dom(5),tech);
            [Wc, ~, ~, K, K2] = ACA(T3, absTol, n(3));
            r(3) = size(K,2);
            IT3 = I;
            JT3 = J;
            
            % Refine n
            breakFlag = 0;
            while ( r(1)*2*sqrt(2) > n(1) )
                n(1) = reffun(n(1));
                breakFlag = 1;
            end
            while ( r(2)*2*sqrt(2) > n(2) )
                n(2) = reffun(n(2));
                for i = 1:3
                    while ( r(i)*2*sqrt(2) > n(i) )
                        n(i) = reffun(n(i));
                    end
                end
                breakFlag = 1;
            end
            while ( r(3)*2*sqrt(2) > n(3) )
                n(3) = reffun(n(3));
                breakFlag = 1;
            end
            
            % Reinitialize after refinement
            if ( breakFlag == 1 )
                happyPhase1 = 0;
                break
                % Proceed to refinement if r < 2
            elseif ( min(r) < 2 )
                break
            end
        end
    end
    
    % Catch the rank zero function:
    if ( size(I,1) == 0 || size(J,1) == 0 || size(K,1) == 0 )
        f.cols = chebfun(zeros([n(1),1]), [dom(1),dom(2)], pref);
        f.rows = chebfun(zeros([n(2),1]), [dom(3),dom(4)], pref);
        f.tubes = chebfun(zeros([n(3),1]), [dom(5),dom(6)], pref);
        f.core = 0;
    else
        %% Refine
        pref.chebfuneps = relTol;
        m = n;
        
        % Check if further refinement is necessary
        Uf = Uc;
        Vf = Vc;
        Wf = Wc;
        [resolvedU, resolvedV, resolvedW] = happinessCheck3(Uf, Vf, Wf, dom, pref, tech);
        
        if ( ~resolvedU )
            if isa(tech,'trigtech')
                m(1) = 2*m(1);
            else
                m(1) = 2*m(1)-1;
            end
        end
        if ( ~resolvedV )
            if isa(tech,'trigtech')
                m(2) = 2*m(2);
            else
                m(2) = 2*m(2)-1;
            end
        end
        if ( ~resolvedW )
            if isa(tech,'trigtech')
                m(3) = 2*m(3);
            else
                m(3) = 2*m(3)-1;
            end
        end
        
        % Add function evaluations and check again
        while ( ~resolvedU || ~resolvedV || ~resolvedW )
            % function handle to evaluate T_f
            ff = @(i,j,k) op(chebX(i,m(1)),chebY(j,m(2)),chebZ(k,m(3)));
            
            % Map the indices from T_c to T_f
            if isa(tech,'trigtech')
                refFactor = m./n;
            else
                refFactor = [0, 0, 0];
                iter = 0;
                while ( min(refFactor) == 0 )
                    iter = iter +1;
                    if ( n(1)*iter-(iter-1) == m(1) )
                        refFactor(1) = iter;
                    end
                    if ( n(2)*iter-(iter-1) == m(2) )
                        refFactor(2) = iter;
                    end
                    if ( n(3)*iter-(iter-1) == m(3) )
                        refFactor(3) = iter;
                    end
                end
            end
            
            if isa(tech,'trigtech')
                ref = @(i, r) r*i-(r-1);
            else
                ref = @(i, r) r*i-(r-1);
            end
            
            % U
            Jr = ref(JT1, refFactor(2));
            Kr = ref(KT1, refFactor(3));
            if ( ~resolvedU )
                Uold = Uf;
                [Ij,Ik] = ind2sub([size(Jr,2),size(Kr,2)],I2);
                if ( vectorize == 0 )
                    nn = [floor(m(1)/2),size(I2,2)];
                    x = zeros([nn(1),1]);
                    x(:,1) = 2:2:m(1);
                    X = repmat(x,1,nn(2));
                    y = zeros([1,nn(2)]);
                    y(1,:) = Jr(Ij);
                    Y = repmat(y,nn(1),1);
                    z = zeros([1,nn(2)]);
                    z(1,:) = Kr(Ik);
                    Z = repmat(z,nn(1),1);
                    Uf(2:2:m(1),1:nn(2)) = ff(X,Y,Z);
                    Uf(1:2:m(1),1:nn(2)) = Uold;
                else
                    Uf = zeros([m(1),size(I2,2)]);
                    for i = 1:m(1)
                        for j = 1:size(I2,2)
                            
                            if mod(i,2) == 0
                                Uf(i,j) = ff(i, Jr(Ij(j)), Kr(Ik(j)));
                            else
                                Uf(i,j) = Uold((i+1)/2,j);
                            end
                        end
                    end
                end
                [resolvedU, ~, ~] = happinessCheck3(Uf, [], [], dom, pref, tech);
                if ( ~resolvedU )
                    if isa(tech,'trigtech')
                        m(1) = 2*m(1);
                    else
                        m(1) = 2*m(1)-1;
                    end
                end
            end
            
            % V
            Ir = ref(IT2, refFactor(1));
            Kr = ref(KT2, refFactor(3));
            if ( ~resolvedV )
                Vold = Vf;
                [Ji,Jk] = ind2sub([size(Ir,2),size(Kr,2)],J2);
                if ( vectorize == 0 )
                    nn = [floor(m(2)/2),size(J2,2)];
                    y = zeros([nn(1),1]);
                    y(:,1) = 2:2:m(2);
                    Y = repmat(y,1,nn(2));
                    x = zeros([1,nn(2)]);
                    x(1,:) = Ir(Ji);
                    X = repmat(x,nn(1),1);
                    z = zeros([1,nn(2)]);
                    z(1,:) = Kr(Jk);
                    Z = repmat(z,nn(1),1);
                    Vf(2:2:m(2),1:nn(2)) = ff(X,Y,Z);
                    Vf(1:2:m(2),1:nn(2)) = Vold;
                else
                    Vf = zeros([m(2),size(J2,2)]);
                    for i = 1:m(2)
                        for j = 1:size(J2,2)
                            if mod(i,2) == 0
                                Vf(i,j) = ff(Ir(Ji(j)), i, Kr(Jk(j)));
                            else
                                Vf(i,j) = Vold((i+1)/2,j);
                            end
                        end
                    end
                end
                [~,resolvedV,~] = happinessCheck3([], Vf, [], dom, pref, tech);
                if ~resolvedV
                    if isa(tech,'trigtech')
                        m(2) = 2*m(2);
                    else
                        m(2) = 2*m(2)-1;
                    end
                end
            end
            
            % W
            Ir = ref(IT3, refFactor(1));
            Jr = ref(JT3, refFactor(2));
            if ( ~resolvedW )
                Wold = Wf;
                [Ki,Kj] = ind2sub([size(Ir,2),size(Jr,2)],K2);
                if ( vectorize == 0 )
                    nn = [floor(m(3)/2),size(K2,2)];
                    z = zeros([nn(1),1]);
                    z(:,1) = 2:2:m(3);
                    Z = repmat(z,1,nn(2));
                    x = zeros([1,nn(2)]);
                    x(1,:) = Ir(Ki);
                    X = repmat(x,nn(1),1);
                    y = zeros([1,nn(2)]);
                    y(1,:) = Jr(Kj);
                    Y = repmat(y,nn(1),1);
                    Wf(2:2:m(3),1:nn(2)) = ff(X,Y,Z);
                    Wf(1:2:m(3),1:nn(2)) = Wold;
                else
                    Wf = zeros([m(3),size(K2,2)]);
                    for i = 1:m(3)
                        for j = 1:size(K2,2)
                            if mod(i,2) == 0
                                Wf(i,j) = ff(Ir(Ki(j)), Jr(Kj(j)), i);
                            else
                                Wf(i,j) = Wold((i+1)/2,j);
                            end
                        end
                    end
                end
                [~,~,resolvedW] = happinessCheck3([], [], Wf, dom, pref, tech);
                if ~resolvedW
                    if isa(tech,'trigtech')
                        m(3) = 2*m(3);
                    else
                        m(3) = 2*m(3)-1;
                    end
                end
            end
        end
        
        
        [~, absTol] = getTol(Uf, pseudoLevel, absTol, dom(2)-dom(1),tech);
        [~, absTol] = getTol(Vf, pseudoLevel, absTol, dom(4)-dom(3),tech);
        [~, absTol] = getTol(Wf, pseudoLevel, absTol, dom(6)-dom(5),tech);
        
        %% Phase 3
        
        % Compute factor matrices
        [QU,~] = qr(Uf,0);
        [I, QUI] = DEIM(QU);
        [QV,~] = qr(Vf,0);
        [J, QVJ] = DEIM(QV);
        [QW,~] = qr(Wf,0);
        [K, QWK] = DEIM(QW);
        
        % Scaling to ensure the factor matrices contain decaying coefficients
        tmpCore  = invtprod(evalTensor(I,J,K,ff,vectorize), QUI,QVJ,QWK);
        colScaling = max(abs(tmpCore),[],[2,3]);
        rowScaling = max(abs(tmpCore),[],[1,3]);
        tubeScaling = squeeze(max(abs(tmpCore),[],[1,2]));
           
        % Store chebfun3 object
        f.cols  = simplify(chebfun(QU*diag(colScaling), [dom(1),dom(2)], pref));
        f.rows  = simplify(chebfun(QV*diag(rowScaling), [dom(3),dom(4)], pref));
        f.tubes = simplify(chebfun(QW*diag(tubeScaling), [dom(5),dom(6)], pref));
        f.cols = simplify(f.cols, pref.chebfuneps, 'globaltol');
        f.rows = simplify(f.rows, pref.chebfuneps, 'globaltol');
        f.tubes = simplify(f.tubes, pref.chebfuneps, 'globaltol');
        f.core  = invtprod(tmpCore,diag(colScaling),diag(rowScaling),diag(tubeScaling));
        
    end
    
    % Sample test
    if ( passSampleTest && restarts <= maxRestarts )
        % Wrap the op with evaluate in case the 'vectorize' flag is on:
        sampleOP = @(x,y,z) evaluate( op, x, y, z, vectorize);
        happy = sampleTest(f, sampleOP, absTol, vectorize);
    else
        happy = 1;
    end
    
    % Restart
    if ( ~happy )
        if ( restarts + 1 == maxRestarts )
            warning('chebfun3f: max number of restarts reached')
            return
        end
        
        % Increase n
        n(1) = floor(sqrt(2)^(floor(2*log2(n(1))) + 1)) + 1;
        n(2) = floor(sqrt(2)^(floor(2*log2(n(2))) + 1)) + 1;
        n(3) = floor(sqrt(2)^(floor(2*log2(n(3))) + 1)) + 1;
        restarts = restarts + 1;
        
        % Ensure r is large enougth for (1,r,r) functions
        if ( r(1) > 1 || r(2) > 1 || r(2) > 1 )
            if r(1) < 3
                r(3) = max(6,2*r(3));
                r(2) = max(6,2*r(2));
            elseif r(2) < 2
                r(1) = max(6,2*r(1));
                r(3) = max(6,2*r(3));
            elseif r(3) < 2
                r(1) = max(6,2*r(1));
                r(2) = max(6,2*r(2));
            end
        end
        
        % Ensure r is large enougth very low-rank functions
        r(1) = max(r(1),3);
        r(2) = max(r(2),3);
        r(3) = max(r(3),3);
    end
end
end

%%
function T = evalTensor(I, J, K, ff,vectorize)
% Evaluate the tensor ff at indices specified by I,J,K

if ( vectorize == 0 ) % we can use the efficient evaluations
    n = [numel(I),numel(J),numel(K)];
    x = zeros([n(1),1,1]);
    x(:,1,1) = I;
    X = repmat(x,1,n(2),n(3));
    y = zeros([1,n(2),1]);
    y(1,:,1) = J;
    Y = repmat(y,n(1),1,n(3));
    z = zeros([1,1,n(3)]);
    z(1,1,:) = K;
    Z = repmat(z,n(1),n(2),1);
    if numel(X) > 0
        T = ff(X,Y,Z);
    else
        T = [];
    end
else % we need for loops as f is not vectorizable
    T = zeros(size(I,2),size(J,2),size(K,2));
    for i = 1:size(I,2)
        for j =1:size(J,2)
            for k = 1:size(K,2)
                T(i,j,k) = ff(I(i),J(j),K(k));
            end
        end
    end
end
end

%%
function [Ac, At, Ar, rowInd, colInd] = ACA(A, tol, maxIter)
% Adaptive Cross Approximation with full pivoting

rowInd = [];
colInd = [];
Aoriginal = A;

for iter = 1:maxIter
    
    [error,I2] = max(abs(A(:)));
    if ( isempty(error) || error < tol )
        Ac = Aoriginal(:,colInd);
        Ar = Aoriginal(rowInd,:)';
        At = Aoriginal(rowInd,colInd);
        return
    end
    
    [I,J] = ind2sub(size(A), I2);
    rowInd = [rowInd, I];
    colInd = [colInd, J];
    
    A = A-A(:,J)*A(I,:)./A(I,J);
end
Ac = Aoriginal(:,colInd);
Ar = Aoriginal(rowInd,:)';
At = Aoriginal(rowInd,colInd);
end

%%
function [indices, UI] = DEIM(U)
% Discrete Empirical Interpolation

indices = [];
[~, I] = max(abs(U(:,1)));
indices = [indices,I];
for l = 2:size(U,2)
    c = U(indices,1:(l-1)) \ U(indices,l);
    r = U(:,l) - U(:,1:(l-1))*c;
    [~, I] = max(abs(r));
    indices = [indices,I];
end

if ( nargout > 1 )
    UI = U(indices,:);
end

end

%%
function ct2 = createCT2(W,data)
% Create temporary chebtech2 object

data.vscale = max(abs(W(:)));
ct2 = chebtech2(W, data);
ct2.coeffs = sum(abs(ct2.coeffs), 2);
end

%%
function [resolvedU, resolvedV, resolvedW] = happinessCheck3(U, V, W, dom, pref, tech)
resolvedU = 1; resolvedV = 1; resolvedW = 1;
% Happiness-check
if numel(U) > 0
    vsclU = max(abs(U(:,1)));
    fiber1Data.hscale = norm(dom(5:6), inf);
    fiber1Data.vscale = vsclU;
    UChebtech = tech.make(U, fiber1Data);
    UChebtech.coeffs = sum(abs(UChebtech.coeffs), 2);
    resolvedU  = happinessCheck(UChebtech, [], ...
        UChebtech.coeffs, [], pref);
end
if numel(V) > 0
    vsclV = max(abs(V(:,1)));
    fiber2Data.hscale = norm(dom(1:2), inf);
    fiber2Data.vscale = vsclV;
    VChebtech = tech.make(V, fiber2Data);
    VChebtech.coeffs = sum(abs(VChebtech.coeffs), 2);
    resolvedV  = happinessCheck(VChebtech, [], ...
        VChebtech.coeffs, [], pref);
end
if numel(W) > 0
    vsclW = max(abs(W(:,1)));
    fiber3Data.hscale = norm(dom(3:4), inf);
    fiber3Data.vscale = vsclW;
    WChebtech = tech.make(W, fiber3Data);
    WChebtech.coeffs = sum(abs(WChebtech.coeffs), 2);
    resolvedW = happinessCheck(WChebtech, [], ...
        WChebtech.coeffs, [], pref);
end

end

%%
function [relTol, absTol] = getTol(M, pseudoLevel, tolOld,domDiff,tech)
% Get suitable tolerances as in Chebfun3 (see
% https://github.com/chebfun/chebfun/issues/1491)

relTol = 2*size(M,1)^(4/5) * pseudoLevel;
vscale = max(abs(M(:)));
cheb = @(i,n) -cos((i-1).*pi/(n-1));
if isa(tech,'trigtech')
    %use equispace points instead
    cheb = @(i,n) -1 +  (i-1)/(n) * 2;
end
points = 1:size(M,1);
points = cheb(points, size(M,1));
gradNorms = zeros([1,size(M,1)]);
for i = 1:size(M,2)
    gradNorms(i) = max(abs(diff(M(:,i)) ./ diff(points)'));
end
gradNorms = max(gradNorms);
absTol = max(max(domDiff.*gradNorms), vscale) * relTol;
absTol = max([absTol, tolOld, pseudoLevel]);
end

%%
function X = initializeIndexRandomly(r, maxVal, restarts)
%  Random initialization of indices by drawing one index in each of r
%  subintervals of equal length

box = floor(maxVal/r);
X = [];
rngprev = rng();
rng(16051821+restarts,'twister');             % pseudorandom
for i = 1:r
    val = i*box + randi(box,1);
    X = [X,val];
end
rng(rngprev);
end

%%
function X = invtprod(X,U,V,W)
% Evaluate X times_1 inv(U) times_2 inv(V) times_3 inv(W) using backslash

n = [size(X,1),size(X,2),size(X,3)];
m = [size(U,1),size(V,1),size(W,1)];
X = reshape(U\reshape(X,[n(1),n(2)*n(3)]),[m(1),n(2),n(3)]);
X = permute(reshape(V\reshape(permute(X,[2,1,3]),[n(2),m(1)*n(3)]),[m(2),m(1),n(3)]),[2,1,3]);
X = permute(reshape(W\reshape(permute(X,[3,2,1]),[n(3),m(2)*m(1)]),[m(3),m(2),m(1)]),[3,2,1]);

end

%%
function vals = evaluate(oper, xx, yy, zz, flag)
% EVALUATE  Wrap the function handle in a FOR loop if the vectorize flag is
% turned on.

if ( flag ==1 )
    if ( isvector(xx) && isvector(yy) && isvector(zz) )
        vals = zeros(size(xx));
        if ( size(xx, 1) == 1 && size(xx, 2) > 1 )
            % Turn rows into columns so that the next for loop works
            % properly.
            xx = xx.';
            yy = yy.';
            zz = zz.';
        end
        for ii = 1: size(xx, 1)
            vals(ii) = oper(xx(ii, 1) , yy(ii, 1), zz(ii, 1));
        end
    else
        vals = zeros(size(xx, 1), size(yy, 2), size(zz, 3));
        for ii = 1:size(xx, 1)
            for jj = 1:size(yy, 2)
                for kk = 1:size(zz, 3)
                    vals(ii, jj, kk) = feval(oper, xx( ii, 1, 1), ...
                        yy(1, jj, 1 ), zz(1, 1, kk));
                end
            end
        end
    end
else % i.e., if (flag == 0)
    vals = feval(oper, xx, yy, zz);  % Tensor or vector of values at cheb3 pts.
    if ( size(vals) ~= size(xx) )
        % Necessary especially when a CHEBFUN2 is made out of a CHEBFUN3,
        % e.g. in CHEBFUN3/STD.
        vals = vals.';
    end
end

end