function f = constructor(f, op, varargin)
%CONSTRUCTOR   CHEBFUN3 constructor.
%   Given a function OP of three variables, this code represents it as a 
%   CHEBFUN3 object. A CHEBFUN3 object is a low-rank representation
%   expressing a function as a trilinear product of a discrete core tensor 
%   and three quasimatrices consisting of univariate functions.
%
%   The algorithm for constructing a CHEBFUN3 has three phases:
%
%   PHASE 1: The first phase attempts to determine the first numerical 
%   separation rank of the function. At the end of this stage we have 
%   candidate 3D pivot locations, 3D pivot values, skeleton slices and 
%   fibers. It also identifies the skeleton columns and rows of each slice.
%   To do so, it uses multivariate adaptive cross approximation (MACA).
%   The output of Phase 1 has the form of block term decomposition (BTD).
%
%   PHASE 2: The second phase attempts to resolve the corresponding 
%   skeleton columns, rows and tubes by sampling along the fibers and 
%   performing MACA on the skeleton. Sampling is increased until the 
%   Chebyshev coefficients of the fibers fall below machine epsilon.
%
%   PHASE 3: The third phase attempts to recompress the BTD form into
%   Tucker format. This step is not an adaptive process. 
%
% See also CHEBFUN2, CHEBFUN3T and CHEBFUN3V.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse the inputs:
[op, dom, pref, vscaleBnd, fiberDim, vectorize, isEqui, fixedRank] = ...
    parseInputs(op, varargin{:});

% Set preferences:
tech            = pref.tech();
tpref           = tech.techPref;
grid            = tpref.minSamples;
maxSample       = tpref.maxLength; % This is the maxSample used in Phase II.
maxSamplePhase1 = 363; % maxSample in Phase I where we explicitly generate 
% order-3 tensors.
prefStruct = pref.cheb3Prefs;
maxRank = prefStruct.maxRank;
pseudoLevel = prefStruct.chebfun3eps;
passSampleTest = prefStruct.sampleTest;

if ( isa(op, 'chebfun3') )     % CHEBFUN3( CHEBFUN3 )
    f = op;
    return
elseif ( isa(op, 'double') )   % CHEBFUN3( DOUBLE )
    f = constructFromDouble(op, dom, pref, isEqui);
    if ( fixedRank )
        % Simplify the rank if requested.
        f = fixTheRank(f , fixedRank);
    end
    return
end

%% Dimension clustering of Bebendorf & Kuske.
% An example where this is important. compare3(@(x,y,z)
% exp(-5*(x+2*y).^2-30*y.^4-7*sin(2*z).^6));
if ( isempty(fiberDim) )
    fiberDim = dimCluster(op, dom, vectorize, pref);
end
if ( fiberDim == 1 )
    op = @(y,z,x) op(x,y,z);
    dom = [dom(3) dom(4) dom(5) dom(6) dom(1) dom(2)];
elseif ( fiberDim == 2 )
    op = @(z,x,y) op(x,y,z);
    dom = [dom(5) dom(6) dom(1) dom(2) dom(3) dom(4)];
end
isHappy = 0;     % We are currently unresolved.
failure = 0;     % Reached max discretization size without being happy.

while ( ~isHappy && ~failure )
    %% Main loop of the constructor
    [xx, yy, zz] = points3D(grid, grid, grid, dom, pref);
    grid2D = grid;
    vals = evaluate(op, xx, yy, zz, vectorize);
    % We have vals(i,j,k) = op(X(i), Y(j), Z(k)), where                 (*)
    %                                           X = chebpts(m,[a,b]);
    %                                           Y = chebpts(n,[c,d]);
    %                                           Z = chebpts(p,[e,g]);
    %                                 [xx, yy, zz] = ndgrid(X, Y, Z);
    % and ndgrid is used here in ``points3D``.
    % If we have used meshgrid, instead of ndgrid, in ``points3D``, then 
    % what holds, instead of (*) above, is 
    % vals(i,j,k) = op(X(j), Y(i), Z(k)) which can be confusing.
    % Users who want to create a chebfun3 from a DISCRETE tensor of values 
    % as input should use ndgrid to generate their input tensor to feed to 
    % chebfun3.
    
    % Does the function blow up or evaluate to NaN?:
    vscale = max(abs(vals(:)));
    if ( isinf(vscale) )
        error('CHEBFUN:CHEBFUN3:constructor:inf', ...
            'Function returned INF when evaluated');
    elseif ( any(isnan(vals(:))) )
        error('CHEBFUN:CHEBFUN3:constructor:nan', ...
            'Function returned NaN when evaluated');
    end
    
    factor = 2*sqrt(2); % Ratio between width of tensor and rank (=no. of pivots).
    % 3D version of CHEBFUN's tolerance: to be used in the ACA for pivot
    % size.
    [relTol, absTol] = getTol3D(xx, yy, zz, vals, grid, dom, pseudoLevel, ...
        vscaleBnd);
    pref.chebfuneps = relTol; % tolerance to be used in happinessCheck.
    
    %% PHASE 1: %%
    % Apply 3D ACA to tensor of values
    [colsValues, rowsValues, pivotVals2D, pivotIndices2D, fibersValues, ...
        pivotVals3D, pivotIndices3D, iFail3D, iFail2D] = completeACA3D(...
        vals, fiberDim, absTol, factor, dom, pref);
    
    strike = 1;
    while ( (iFail3D || iFail2D) && grid < maxSamplePhase1 && strike < 3)
        % Refine sampling on tensor grid:
        if ( iFail3D )
            grid = gridRefinePhase1(grid, pref);
            [xx, yy, zz] = points3D(grid, grid, grid, dom, pref);
        elseif ( iFail2D )
            grid2D = gridRefinePhase1(grid2D, pref);
            [xx, yy, zz] = points3D(grid2D, grid2D, grid, dom, pref);
        end
        
        vals = evaluate(op, xx, yy, zz, vectorize); % Resample
        vscale = max(abs(vals(:)));
        % New tolerance:
        [relTol, absTol] = getTol3D(xx, yy, zz, vals, grid, dom, ...
            pseudoLevel, vscaleBnd);
        pref.chebfuneps = relTol;
        
        % New 3D ACA:
        [colsValues, rowsValues, pivotVals2D, pivotIndices2D, ...
            fibersValues, pivotVals3D, pivotIndices3D, iFail3D, iFail2D] ...
            = completeACA3D(vals, fiberDim, absTol, factor, dom, pref);
        % If the function is 0+noise then stop after three strikes.
        if ( abs(pivotVals3D(1))<1e4*absTol/vscale )
            strike = strike + 1;
        end
    end
    
    % If the rank of the function is above maxRank then stop.
    if ( numel(pivotVals3D) > maxRank )
        warning('CHEBFUN:CHEBFUN3:constructor:rank', ...
            'Not a low-rank function.');
        failure = 1;
    end
    
    % Check if the pivot rows, columns and tubs are resolved.
    [resolvedTubs, resolvedCols, resolvedRows] = ...
        happinessCheckBTD(fibersValues, colsValues, rowsValues, dom, ...
        pref, tech);
    isHappy = resolvedCols & resolvedRows & resolvedTubs;
    
    %% Prepare for Phase 2:
    [gridX, gridY, gridZ] = size(vals);
    sepRank = numel(pivotVals3D); % first separation rank
    
    % Check whether the function is zero:
    if ( length(pivotVals3D) == 1 && pivotVals3D == 0 )
        isHappy = 1;
    else
        % Find the location of 3D and 2D pivot values:
        pivPos3D = [ xx(pivotIndices3D(:, 1), 1, 1), ...
            yy(1, pivotIndices3D(:, 2), 1).', ...
            squeeze(zz(1, 1, pivotIndices3D(:, 3)))];
        PI3D = pivotIndices3D; % This will be overwritten in Phase 2, if 
        % necessary.
        for k = 1:sepRank
            pivPos2D{k} = [ xx(pivotIndices2D{k}(:, 1), 1, 1), ...
                yy(1,pivotIndices2D{k}(:, 2), 1).'];
        end
        PI2D = pivotIndices2D;
        for k=1:sepRank
            diagValues2D{k} = diag(1./pivotVals2D{k});
        end
    end
    %% PHASE 2:
    % Separation rank of f is already detected and the skeleton columns, 
    % rows and tubes are also found in Phase 1, i.e., a BTD is computed. 
    % However, if f is not yet resolved, we refine our sampling of f only 
    % along the skeleton tubes and skeleton columns and rows of the 
    % skeleton slices.
    grid_Phase1 = grid;
    failure = 0;
    while ( (~isHappy) && (~failure) )
        m = gridX; 
        n = gridY;
        p = gridZ;
        % Identify direction/directions that need refinement:
        if ( ~resolvedTubs )
            % Double sampling along skeleton tubes
            [p, nesting] = gridRefinePhase2(p, pref);
            Z = mypoints(p, dom(5:6), pref);
            fibersValues = zeros(p, sepRank);
            for k=1:sepRank
                [xx, yy, zz] = ndgrid(pivPos3D(k, 1), pivPos3D(k, 2), Z);
                fibersValues(:,k) = squeeze(evaluate(op, xx, yy, zz, vectorize));
                % An alternative is the following, which does
                % not form xx, yy and zz, but might be easier to read.
		% fibersValues(:,k) = evaluate(op, repmat(pivPos3D(k, 1),p,1),
		% repmat(pivPos3D(k, 2),p,1), Z, vectorize ).';
            end
            
            % Find location of pivots on new grid  (using nesting property).
            PI3D(:, 3) = nesting(PI3D(:, 3)); % The first column of PP contains 
            % z indices and should therefore be updated.
        else
            fibersValues = zeros(p, sepRank);
            Z = mypoints(p, dom(5:6), pref);
            for k=1:sepRank
                [xx, yy, zz] = ndgrid(pivPos3D(k, 1), pivPos3D(k, 2), Z);
                fibersValues(:, k) = squeeze(evaluate(op, xx, yy, zz, vectorize));
            end
        end
        
        if ( ~resolvedCols && resolvedRows )
            % Double sampling along the column slices only, i.e., in x. 
            [m, nesting] = gridRefinePhase2(m, pref); 
            X = mypoints(m, dom(1:2), pref);
            Y = mypoints(n, dom(3:4), pref);
            colsValues = {};
            for k=1:sepRank
                [xx, yy, zz] = ndgrid(X, pivPos2D{k}(:,2), pivPos3D(k, 3));
                colsValues{k} = evaluate(op, xx, yy, zz, vectorize);
            end
            
            rowsValues = {};
            for k=1:sepRank
                [xx, yy, zz] = ndgrid(pivPos2D{k}(:,1), Y, pivPos3D(k, 3));
                rowsValues{k} = evaluate(op, xx, yy, zz, vectorize).';
            end
            
            % Find location of 3D pivots on new grid  (using nesting property).
            PI3D(:, 1) = nesting(PI3D(:, 1)); % The 1st column of PI3D contains x
            % indices and should therefore be updated.
	    % Find location of 2D pivots on new grid  (using nesting property).
            for k = 1:sepRank
                PI2D{k}(:,1) = nesting(PI2D{k}(:, 1));
            end
            
        elseif ( resolvedCols && ~resolvedRows )
            [n, nesting] = gridRefinePhase2(n, pref);
            X = mypoints(m, dom(1:2), pref);
            Y = mypoints(n, dom(3:4), pref);

            colsValues = {};
            for k=1:sepRank
                [xx, yy, zz] = ndgrid(X, pivPos2D{k}(:, 2), pivPos3D(k, 3));
                colsValues{k} = evaluate(op, xx, yy, zz, vectorize);
            end
            
            rowsValues = {};
            for k=1:sepRank
                [xx, yy, zz] = ndgrid(pivPos2D{k}(:, 1), Y, pivPos3D(k, 3));
                rowsValues{k} = evaluate(op, xx, yy, zz, vectorize).';
            end
            
            PI3D(:, 2) = nesting(PI3D(:, 2)); % The 3rd column of PP 
            % contains z indices and should therefore be updated.
            for k = 1:sepRank
                PI2D{k}(:, 2) = nesting(PI2D{k}(:, 2));
            end
            
        elseif ( ~resolvedCols && ~resolvedRows )
            [m, nesting1] = gridRefinePhase2(m, pref);
            [n, nesting2] = gridRefinePhase2(n, pref);
            X = mypoints(m, dom(1:2), pref);
            Y = mypoints(n, dom(3:4), pref);
            
            colsValues = {};
            for k=1:sepRank
                [xx, yy, zz] = ndgrid(X, pivPos2D{k}(:, 2), pivPos3D(k, 3));
                colsValues{k} = evaluate(op, xx, yy, zz, vectorize);
            end
            
            rowsValues = {};
            for k=1:sepRank
                [xx, yy, zz] = ndgrid( pivPos2D{k}(:, 1), Y, pivPos3D(k, 3));
                rowsValues{k} = evaluate(op, xx, yy, zz, vectorize).';
            end
            
            PI3D(:, 1) = nesting1(PI3D(:, 1));
            PI3D(:, 2) = nesting2(PI3D(:, 2));
            for k = 1:sepRank
                PI2D{k}(:,1) = nesting1(PI2D{k}(:, 1));
                PI2D{k}(:,2) = nesting2(PI2D{k}(:, 2));
            end
            
        end
        
        % Now, apply ACA to the skeletons only. This is analogous to Phase 
        % 2 of Chebfun2. N.B.: Two levels of ACA on skeletons are needed. 
        % First a 3D ACA should update cols{kk}, rows{kk} and tubes(:, kk). 
        % Then, a 2D ACA should be applied to cols{kk} and rows{kk} so that 
        % they are finalized and ready for use in updating next cols, rows, and
        % tubes.
        for kk = 1:sepRank            
            % N.B.: This for-loop is designed so that at each step, the 3D
            % ACA is applied and then the necessary 2D ACA for each cols 
            % and rows is applied in the NEXT iteration. 
            % The 1st slice of BTD needs no 3D ACA, but needs a 2D ACA and
            % this is why the 2D ACA appears earlier than the 3D ACA in the
            % following lines. On the other hand, at the last iteration, 
            % when kk = sepRank, only 2D ACA updates cols{kk} and rows{kk}.
            % The 3D ACA needed on cols{sepRank} and rows{sepRank} is 
            % actually done in the previous iteration, i.e., when kk has
            % been equal to sepRank - 1. So, no further 3D ACA should affect 
            % cols{sepRank} and rows{sepRank} in that last step.
            
            jj=kk;
            nn2D = numel(pivotVals2D{jj});
            if ( ~resolvedCols || ~resolvedRows )            
            for iii = 1:nn2D-1                                
            colsValues{jj}(:, iii+1:end) = colsValues{jj}(:, iii+1:end) -...
                (colsValues{jj}(:, iii)*diagValues2D{jj}(iii,iii)) * ...
                (rowsValues{jj}(PI2D{jj}(iii+1:nn2D,2), iii )).';
                
            rowsValues{jj}(:, iii+1:end) = rowsValues{jj}(:, iii+1:end) -...
                ((colsValues{jj}(PI2D{jj}(iii+1:nn2D,1), iii ) * ...
                diagValues2D{jj}(iii,iii))*(rowsValues{jj}(:, iii).')).';
            end
            end

            % Apply 3D ACA just on the the skeleton columns, rows and tubes of 
            % the BTD to update them. This is a 3D analogue of Phase 2 of 
            % chebfun2.
            
            % Step 1: update colsValues
            if ( ~resolvedCols || ~resolvedRows )
                for ii = kk+1:sepRank
                    colsValues{ii} = colsValues{ii} - (colsValues{kk} * ...
                        diagValues2D{kk})*rowsValues{kk}(PI2D{ii}(:,2),:).' ...
                        *(fibersValues(PI3D(ii, 3),kk) ./ pivotVals3D(kk));
                end
                
                % Step 2: update rowsValues
                for ii = kk+1:sepRank
                    rowsValues{ii} = rowsValues{ii} - ...
                        (colsValues{kk}(PI2D{ii}(:,1), :) * ...
                        diagValues2D{kk} * rowsValues{kk}.').' * ...
                        (fibersValues(PI3D(ii, 3), kk) ./ pivotVals3D(kk));
                end
            end
            
            % Step 3: update fibersValues
            for ii = kk+1:sepRank
                fibersValues(:, ii) = fibersValues(:, ii) - ...
                    (colsValues{kk}(PI3D(ii, 1), :) * diagValues2D{kk} ...
                    * rowsValues{kk}(PI3D(ii,2),:).') ...
                    * (fibersValues(:,kk) ./ pivotVals3D(kk));
            end
        end
        
        %% 3D HAPPINESSCHECK again only on non-happy fibers:
        fiber1Data.hscale = norm(dom(1:2), inf);
        fiber1Data.vscale = max(abs(colsValues{1}(:)));
        if ( ~resolvedCols )
            colsChebtech = tech.make(sum(colsValues{1}, 2), fiber1Data);
            resolvedCols  = happinessCheck(colsChebtech ,[], ...
                sum(colsValues{1}, 2), [], pref);
        end

        fiber2Data.hscale = norm(dom(3:4), inf);
        fiber2Data.vscale = max(abs(rowsValues{1}(:)));
        if ( ~resolvedRows )
            rowsChebtech = tech.make(sum(rowsValues{1}, 2), fiber2Data);
            resolvedRows  = happinessCheck(rowsChebtech, [], ...
                sum(rowsValues{1}, 2), [], pref);
        end

        fiber3Data.hscale = norm(dom(5:6), inf);
        fiber3Data.vscale = max(abs(fibersValues(:)));
        if ( ~resolvedTubs )
            fiber3Chebtech = tech.make(sum(fibersValues, 2), fiber3Data);
            resolvedTubs  = happinessCheck(fiber3Chebtech, [], ...
                sum(fibersValues, 2), [], pref);
        end        
        isHappy = resolvedCols & resolvedRows & resolvedTubs;
        
        if ( ~resolvedCols )
            gridX = gridRefinePhase2(gridX, pref);
        else
            gridX = m;    
        end
        if ( ~resolvedRows )
            gridY = gridRefinePhase2(gridY, pref);
        else
            gridY = n;
        end
        if ( ~resolvedTubs )
            gridZ = gridRefinePhase2(gridZ, pref);
        else
            gridZ = p;
        end
        
        % STOP if degree is over maxLength:
        if ( max([gridX, gridY, gridZ]) >= maxSample )
            warning('CHEBFUN:CHEBFUN3:constructor:notResolved', ...
                'Unresolved with maximum CHEBFUN3 length: %u.', maxSample);
            failure = 1;
        end
        grid = min([m, n, p]);
    end
    %% End of Phase II
    % If f is numerically zero, artificially set the columns, rows and 
    % fibers to zero.
    if ( (norm(colsValues{1}) == 0) || (norm(rowsValues{1}) == 0) || ...
            (norm(fibersValues) == 0) )
        fibersValues = zeros(grid, 1);
        colsValues{1} = zeros(grid, 1);
        rowsValues{1} = zeros(grid, 1);
        diagValues2D{1} = 0;        
        pivotVals3D = Inf;
        isHappy = 1;
    end
    
    %% PHASE 3:
    % BTD2Tucker compression. N.B.: This is a non-adaptive step in the
    % sense that it works with a fixed BTD.
    if ( (isHappy) || (failure) )
        dom2D = dom(1:4); % The first 4 entries correspond to slices
        [core, colsValues, rowsValues] = btd2tucker(colsValues, ...
            rowsValues, diagValues2D, pivotVals3D, dom2D, pref, absTol);
        % Developer Note: To check the validity of the new representation 
        % try the following:
	% TuckerApprox = txm(txm(txm(core, colsValues, 1), rowsValues, 2),
	% fibersValues, 3);
        % norm123 = norm( vals(:) -  TuckerApprox(:) )
    end
    
    % Construct a CHEBFUN3 object, simplify skeletons and call a
    % sampleTest.
    f.cols = chebfun(colsValues, dom(1:2), pref);
    f.rows = chebfun(rowsValues, dom(3:4), pref);
    f.tubes = chebfun(fibersValues, dom(5:6), pref);    
    f.cols = simplify(f.cols, pref.chebfuneps, 'globaltol');
    f.rows = simplify(f.rows, pref.chebfuneps, 'globaltol');
    f.tubes = simplify(f.tubes, pref.chebfuneps, 'globaltol');
    
    f.core = core;
    f.domain = dom;
    
    % Sample Test
    if ( passSampleTest )
        % Wrap the op with evaluate in case the 'vectorize' flag is on: 
        sampleOP = @(x,y,z) evaluate( op, x, y, z, vectorize);
        pass = sampleTest(f, sampleOP, absTol, vectorize);
        if ( ~pass )
            isHappy = 0;
        end
    end
    
    if ( ~isHappy )
        grid = gridRefinePhase1(grid_Phase1, pref);
    end
end

% Permute back the object if necessary:
if ( fiberDim == 1 )
    f = permute(f, [3 1 2]);
    % [3 1 2] brings the original [2 3 1] permutation back to [1 2 3].
elseif ( fiberDim == 2 )
    f = permute(f, [2 3 1]);
    % [2 3 1] brings the original [3 1 2] permutation back to [1 2 3].
end

if ( all(fixedRank) )
    % Truncate if requested.
    f = fixTheRank(f , fixedRank);
end

end
%% End of constructor

function [core, colsTucker, rowsTucker] = btd2tucker(colsValues, ...
    rowsValues, diagValues2D, pivotVals3D, dom2, pref, absTol)
allCols = []; 
allRows = []; 
allDiags = [];
nn = numel(pivotVals3D);
sizeIndex = zeros(nn, 1);
pseudoLevel = pref.cheb3Prefs.chebfun3eps;

globalTolCols = absTol; % Like the tolerance used to stop adding new fibers
globalTolRows = absTol; % in Phase 1 of 3D ACA.

for kkk = 1:nn
    sizeIndex(kkk) = size(rowsValues{kkk}, 2);
    allCols = [allCols, colsValues{kkk}];
    allRows = [allRows, rowsValues{kkk}];
    allDiags = [allDiags; diag(diagValues2D{kkk})];
end
sizeIndex = cumsum([0; sizeIndex]);
        
% Compress U and V:
if size(allCols,2)>1
    [colsTucker, Su, Vu] = chebfun2ACA(allCols, dom2, pref, ...
        [], globalTolCols); % the 0-flag is used because the tolerance in 
    % chebfun2ACA should be computed differently from the case where the 
    % input matrix is actually vals in chebpts. Here, columns of allCols do
    % not depend on chebpts in that same way because of the way we created 
    % allCols is just by appending all the columns of colsValues. The same 
    % is true for allRows. If we compute tol in the same way as function 
    % values, then the accuracy unreasonably decreases compared to the 
    % slice decomposition.
    Vu = Vu*diag(1./Su);
else
    colsTucker = allCols; Vu = 1;
end

if size(allRows,2)>1
    [rowsTucker,Sv,Vv] = chebfun2ACA(allRows, dom2, pref, ...
        [], globalTolRows);
    Vv = Vv*diag(1./Sv);
else
    rowsTucker = allRows; 
    Vv = 1;
end

% Form the core tensor:
core = zeros(size(colsTucker, 2), size(rowsTucker, 2), nn);
for kkk = 1:nn
    core(:,:,kkk) = Vu(sizeIndex(kkk)+1:sizeIndex(kkk+1), :).' * ...
        diag(allDiags(sizeIndex(kkk)+1:sizeIndex(kkk+1))) * ...
        Vv(sizeIndex(kkk)+1:sizeIndex(kkk+1), :)./pivotVals3D(kkk);
end

if ( (max(abs(allCols(:))) == 0) || (max(abs(allRows(:))) == 0) ) 
    % zero input
    core = 0;
end

end

%%
function tol = GetTol2D(xx, yy, vals, dom, pseudoLevel)
% GETTOL2D     Calculate a tolerance like the one in the Chebfun2 constructor.
%
%  This is the 2D analogue of the tolerance employed in the chebtech
%  constructors. It is based on a finite difference approximation to the
%  gradient, the size of the approximation domain, the internal working
%  tolerance, and an arbitrary (2/3) exponent. 

[m, n] = size(vals);
grid = max(m, n);

if ( n == 1 ) % f is probably rank-1, i.e., vals depend only on one of the 
              % variables.
    df = diff(vals, 1, 1) ./ diff(xx, 1, 1); 
    Jac_norm = max( abs(df(:)) );
else % f probably has a higher rank than 1
    % Remove some edge values so that df_dx and df_dy have the same size. 
    % xx is generated by ndgrid, i.e., xx changes in the first mode:
    dfdx = diff(vals(:,1:n-1), 1, 1) ./ diff(xx(:, 1:n-1), 1, 1); 
    % yy is generated by ndgrid, i.e., yy changes row-wise (2nd mode):
    dfdy = diff(vals(1:m-1,:), 1, 2) ./ diff(yy(1:m-1, :), 1, 2);
    % An approximation for the norm of the gradient over the whole domain.
    Jac_norm = max( max( abs(dfdx(:)), abs(dfdy(:)) ) );
end
vscale = max(abs(vals(:)));
tol = grid^(4/5) * max(abs(dom(:))) * max(Jac_norm^(5/8), vscale) * pseudoLevel;

end
%%
function tol = GetTol2Dv2(xx, yy, vals, dom, pseudoLevel)
% GETTOL2D   Calculate a tolerance for the Chebfun2 constructor.
%
%  This is the 2D analogue of the tolerance employed in the chebtech
%  constructors. It is based on a finite difference approximation to the
%  gradient, the size of the approximation domain, the internal working
%  tolerance, and an arbitrary (4/5) exponent. 

[m, n] = size(vals); 
grid = max(m, n);

% Remove some edge values so that df_dx and df_dy have the same size. 
% xx is generated by ndgrid, i.e., xx changes in the first mode:
dfdx = diff(vals(:, 1:n-1), 1, 1) ./ diff(xx(:, 1:n-1), 1, 1);
% yy is generated by ndgrid, i.e., yy changes row-wise (2nd mode):
dfdy = diff(vals(1:m-1, :), 1, 2) ./ diff(yy(1:m-1, :), 1, 2);

% An approximation for the norm of the gradient over the whole domain.
Jac_norm = max(max(abs(dfdx(:)), abs(dfdy(:))));
vscale = max(abs(vals(:)));
tol = grid^(4/5) * max(abs(dom(:))) * max(Jac_norm, vscale) * pseudoLevel;
end

%%
function [col, pivotVals, row, pivotLoc] = chebfun2ACA(op, dom, pref, ...
    flag_funVals, globalTol)

pseudoLevel = pref.cheb3Prefs.chebfun3eps;
if ( ~isempty(globalTol) )
    tol = globalTol;
else
    if flag_funVals
        [xx2D, yy2D] = points2D(size(op,1), size(op,2), dom, pref); 
        % ndgrid is used here.
        tol = GetTol2D(xx2D, yy2D, op, dom, pseudoLevel);
    else
        vscale = max(abs(op(:)));
        tol = length(op)^(3/4) * max(abs(dom(:))) * vscale * pseudoLevel;
    end
end

% Perform GE with complete pivoting:
[pivotVals, pivotLoc, row, col] = completeACA2D(op, tol, 0); % factor = 0, 
% because we want the ACA to be applied even if op is not low-rank.
% In contrast to Chebfun2, we now have op = col*diag(1./pivotVals)*row'.
end

%%
function [col, pivotVals, row, pivotLoc, ifail2D] = chebfun2ACAv2(op, ...
    tol, factor)
% Perform GE with complete pivoting:

if factor ~= 0
    % FACTOR in the 3D steps is either 0 (in case of constructionFromDoubles) 
    % or 2sqrt(2) otherwise. For 2D steps however, we are happy with 
    % FACTOR = 0 or 2 as in Chebfun2. This IF conditional, makes it
    % possible to rewrite FACTOR in 2D steps, but at the same time keeping 
    % it zero fro constructionFromDoubles.
    factor = 2;
end
[pivotVals, pivotLoc, row, col, ifail2D] = completeACA2D(op, tol, factor); 
end

%%
function [pivotValue, pivotElement, rows, cols, ifail2D] = ...
    completeACA2D(A, tol, factor) 
% 2D ACA with complete pivoting which is the continuous analogue of 
% Gaussian elimination with complete pivoting.
% We attempt to adaptively find the numerical rank of function in the 2D
% level. This is _almost_ the same as the one in chebfun2/constructor.


% Set up output variables.
[nx, ny] = size(A);
width = min(nx, ny);        % Use to tell us how many pivots we can take.
pivotValue = zeros(1);      % Store an unknown number of Pivot values.
pivotElement = zeros(1, 2); % Store (j,k) entries of pivot location.
ifail2D = 1;                  % Assume we fail.

% Main algorithm
zRows = 0;                  % count number of zero cols/rows.
[infNorm, ind] = max(abs(reshape(A, numel(A), 1)));
[row, col] = chebfun3.myind2sub(size(A) , ind);

% Bias toward diagonal for square matrices (see reasoning below):
if ( ( nx == ny ) && ( max(abs(diag(A))) - infNorm) > -tol )
    [infNorm, ind] = max(abs(diag(A)));
    row = ind;
    col = ind;
end

scl = infNorm;
% If the function is the zero function
if ( scl == 0 )
    pivotValue = 0;
    cols = 0;
    rows = 0;
    ifail2D = 0;
else
    cols(:,1) = zeros(size(A, 1), 1);
    rows(1,:) = zeros(1, size(A, 2));
end

while ( ( infNorm > tol ) && ( zRows < width / factor) ...
        && ( zRows < min(nx, ny) ) )

    cols(:, zRows+1) = A(:, col);             % Extract skeleton columns
    rows(zRows+1, :) = A(row, :);             % Extract skeleton rows
    PivVal = A(row, col);
    A = A - cols(:, zRows+1)*(rows(zRows+1,:)./PivVal); % One step of GE
    
    % Keep track of progress.
    zRows = zRows + 1;                       % One more row is interpolated
    pivotValue(zRows) = PivVal;              % Store value of 2D pivot
    pivotElement(zRows, :)=[row col];        % Store index of 2D pivot
    
    % Find value and index of next 2D pivot
    [infNorm, ind] = max(abs(A(:))); % Slightly faster
    [row, col] = chebfun3.myind2sub(size(A), ind);
    
    % Have a bias towards the diagonal of A, so that it can be used as a test
    % for nonnegative definite functions. (Complete GE and Cholesky are the
    % same as nonnegative definite functions have an absolute maximum on the
    % diagonal, except there is the possibility of a tie with an off-diagonal
    % absolute maximum. Bias toward diagonal maxima to prevent this.)
    if ( ( nx == ny ) && ( max(abs(diag(A))) - infNorm) > -tol )
        [infNorm, ind] = max(abs(diag(A)));
        row = ind;
        col = ind;
    end
end

if ( infNorm <= tol )
    ifail2D = 0;                               % We didn't fail in 2D ACA
end
if ( zRows >= (width/factor) )
    ifail2D = 1;                               % We did fail in 2D ACA
end

rows = rows.';                               % To unify all the columns, 
                                             % rows and tubes, store 
                                             % skeleton rows also as column 
                                             % vectors.
end

%%
function [relTol, absTol] = getTol3D(xx, yy, zz, vals, grid, dom, ...
    pseudoLevel, vscaleBnd)
%See https://github.com/chebfun/chebfun/issues/1491

relTol = 2*grid^(4/5) * pseudoLevel; % this should be vscale and hscale invariant
vscale = max(abs(vals(:)));

if ( isempty(vscaleBnd) )
    [m,n,p] = size(vals);
    % Remove some edge values so that df_dx, df_dy and df_dz have the same size. 
    % xx changes in the first mode:
    df_dx = diff(vals(:, 1:n-1, 1:p-1), 1, 1) ./ diff(xx(:, 1:n-1, 1:p-1), 1, 1);
    % yy changes row-wise (2nd mode):
    df_dy = diff(vals(1:m-1, :, 1:p-1), 1, 2) ./ diff(yy(1:m-1, :, 1:p-1), 1, 2);
    % zz changes tube-wise (3rd mode):
    df_dz = diff(vals(1:m-1, 1:n-1, :), 1, 3) ./ diff(zz(1:m-1, 1:n-1, :), 1, 3);
    J = max(max(abs(df_dx),abs(df_dy)), abs(df_dz));
    Jac_norm = max(J(:)); % An approximation for the norm of the Jacobian over
                          % the whole domain.
    absTol =  max(abs(dom(:))) * max(Jac_norm, vscale) * relTol;
else % If we have a binary operation and vscale of both operand are 
    % available in the vector vscaleBnd. 
    kappa = sum(vscaleBnd)/vscale; % condition number of the plus operation 
    % f + g is ||f|| + ||g|| / ||f+g||. 
    absTol = max(abs(dom(:))) * max(kappa, vscale) * relTol;
    % Since kappa is big, if g ~= -f, this helps us find out cases prone to 
    % cancellation errors, i.e., in this case we expect less from the
    % constructor. How less is exactly according to how much cancellation
    % there might be. Check e.g. the value of kappa for f+g with
    % f = chebfun3(@(x,y,z) x ); g = f-0.9*f;
    % and
    % f = chebfun3(@(x,y,z) x ); g = f-0.99999999*f;
end

end

%%
function x = mypoints(n, dom, pref)
% Get the sample points that correspond to the right grid for a particular
% technology.

% What tech am I based on?:
tech = pref.tech();

if ( isa(tech, 'chebtech2'))
    x = chebpts(n, dom, 2);
elseif ( isa(tech, 'chebtech1'))
    x = chebpts(n, dom, 1);
elseif ( isa(tech, 'trigtech'))
    x = trigpts(n, dom);
else
    error('CHEBFUN:CHEBFUN3:constructor:mypoints:techType', ...
        'Unrecognized technology');
end

end

%%
function [xx, yy] = points2D(m, n, dom, pref)
% Get the sample points that correspond to the right grid for a particular
% technology.

% What tech am I based on?:
tech = pref.tech();

if ( isa(tech, 'chebtech2') )
    x = chebpts(m, dom(1:2), 2);   % x grid.
    y = chebpts(n, dom(3:4), 2);   % y grid.
    [xx, yy] = ndgrid(x, y);
elseif ( isa(tech, 'chebtech1'))
    x = chebpts(m, dom(1:2), 1);   % x grid.
    y = chebpts(n, dom(3:4), 1);   % y grid.
    [xx, yy] = ndgrid(x, y);
elseif ( isa(tech, 'trigtech'))
    x = trigpts(m, dom(1:2));      % x grid.
    y = trigpts(n, dom(3:4));      % y grid.
    [xx, yy] = ndgrid(x, y);
else
    error('CHEBFUN:CHEBFUN3:constructor:points2D:tecType', ...
        'Unrecognized technology');
end

end

%%
function [xx, yy, zz] = points3D(m, n, p, dom, pref)
% Get the sample points that correspond to the right grid for a particular
% technology.

% What tech am I based on?:
tech = pref.tech();

if ( isa(tech, 'chebtech2') )
    x = chebpts(m, dom(1:2), 2);   % x grid.
    y = chebpts(n, dom(3:4), 2);   % y grid.
    z = chebpts(p, dom(5:6), 2);   % z grid.
    [xx, yy, zz] = ndgrid(x, y, z); 
elseif ( isa(tech, 'chebtech1') )
    x = chebpts(m, dom(1:2), 1);   % x grid.
    y = chebpts(n, dom(3:4), 1);   % y grid.
    z = chebpts(p, dom(5:6), 1);   % z grid.
    [xx, yy, zz] = ndgrid(x, y, z); 
elseif ( isa(tech, 'trigtech') )
    x = trigpts(m, dom(1:2));      % x grid.
    y = trigpts(n, dom(3:4));      % y grid.
    z = trigpts(p, dom(5:6));      % z grid.
    [xx, yy, zz] = ndgrid(x, y, z);
else
    error('CHEBFUN:CHEBFUN3:constructor:points3D:tecType', ...
        'Unrecognized technology');
end
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

%%
function grid = gridRefinePhase1(grid, pref)
% Grid refinement strategy in Phase 1: It increases sampling by a factor 
% of sqrt(2) on tensor grid for all TECHS.

% What tech am I based on?:
tech = pref.tech();

% What is the next grid size?
if ( isa(tech, 'chebtech2') )
    % Increase sampling by a factor of sqrt(2) on tensor grid:
    grid = floor(sqrt(2)^(floor(2*log2(grid)) + 1)) + 1;
elseif ( isa(tech, 'trigtech') )
    % Increase sampling to an even factor of sqrt(2) on tensor grid:
    grid = round((floor(sqrt(2)^( floor(2*log2(grid)) + 1)) + 1)/2)*2;
    %N.B.: X = round(X./n).*n rounds to a multiple of n, in this case n = 2
elseif ( isa(tech, 'chebtech1') )
    grid = floor(sqrt(2)^(floor(2*log2(grid)) + 1)) + 1;
else
    error('CHEBFUN:CHEBFUN3:constructor:gridRefinePhase1:techType', ...
        'Technology is unrecognized.');
end

end

%%
function [grid, nesting] = gridRefinePhase2(grid, pref)
% Grid refinement strategy for tech in Phase 2

% What tech am I based on?:
tech = pref.tech();

% What is the next grid size?
if ( isa(tech, 'chebtech2') )
    % Double sampling on tensor grid:
    grid = 2*grid-1;
    nesting = 1:2:grid;
elseif ( isa(tech, 'trigtech') )
    % Double sampling on tensor grid:
    grid = 2*grid;
    nesting = 1:2:grid;
elseif ( isa(tech, 'chebtech1' ) )
    % Triple sampling on tensor grid:
    grid = 3 * grid; 
    nesting = 2:3:grid; 
else
    error('CHEBFUN:CHEBFUN3:constructor:gridRefinePhase2:techType', ...
        'Technology is unrecognized.');
end

end

%%
function fiberDim = dimCluster(op, dom, vectorize, pref)
% Choose the right dimension for clustering in Phase 1. The 3 steps are as
% follows: 
% 1) Sample OP at a small tensor.
% 2) Compare all the three modal ranks and find which rank might be the
%    smallest. We are trying to avoid using lots of points that we need in 
%    the technique by Bebendorf and Kuske.
% 3) If there are more than one minimal rank, find the variable that needs
%    more coefficients to resolve.

grid = 10;
[xx, yy, zz] = points3D(grid, grid, grid, dom, pref);
vals = evaluate(op, xx, yy, zz, vectorize);

% Method 1: Using SVD and 1D Chebfun:
F1 = chebfun3.unfold(vals, 1);
F2 = chebfun3.unfold(vals, 2);
F3 = chebfun3.unfold(vals, 3);
rX = rank(F1.'); % Transpose to have longer columns
rY = rank(F2.'); % and therefore being faster
rZ = rank(F3.'); % in MATLAB.
r = [rX, rY, rZ];
[ignored, ind] = find(r==min(r));
if ( numel(ind) == 1 )
    fiberDim = ind;
    return
end
% More than one minimum ranks exist:
tol = 1e-3;
lenX = length(simplify(chebfun(F1(:, 1)), tol));
lenY = length(simplify(chebfun(F2(:, 1)), tol));
lenZ = length(simplify(chebfun(F3(:, 1)), tol));
len = [lenX, lenY, lenZ];
if ( numel(ind) == 2 )
    [ignored, index] = max([len(ind(1)), len(ind(2))]);
    fiberDim = ind(index);
else % numel(ind) = 3
    [ignored, fiberDim] = max(len);
end

%%
% % Method 2: Using Chebfun2:
% F1 = chebfun2(chebfun3.unfold(vals, 1), 
%               [dom(1:2) min(dom(3), dom(5)) max(dom(4), dom(6))]);
% F2 = chebfun2(chebfun3.unfold(vals, 2), 
%               [dom(3:4) min(dom(1), dom(5)) max(dom(2), dom(6))]);
% F3 = chebfun2(chebfun3.unfold(vals, 3), 
%               [dom(5:6) min(dom(1), dom(3)) max(dom(2), dom(4))]);
% rX = numel(F1.pivotValues);
% rY = numel(F2.pivotValues);
% rZ = numel(F3.pivotValues);
% r = [rX, rY, rZ];
% [ignored, ind] = find(r==min(r));
% if ( numel(ind) == 1 )
%     fiberDim = ind;
%     return
% end
% % More than one minimum ranks exist. So, find the variable that needs
% % largest number of coefficients.
% F1 = simplify(F1, 1e-3);
% F2 = simplify(F2, 1e-3);
% F3 = simplify(F3, 1e-3);
% [ignored, lenX] = length(F1); % (simplified) degree needed in x
% [ignored, lenY] = length(F2); % (simplified) degree needed in y
% [ignored, lenZ] = length(F3); % (simplified) degree needed in z
% len = [lenX, lenY, lenZ];
% if ( numel(ind) == 2 )
%     [ignored, index] = max([len(ind(1)), len(ind(2))]);
%     fiberDim = ind(index);
% else % numel(ind) = 3
%     [ignored, fiberDim] = max(len);
% end
end

%%
function [resolvedFibers, resolvedSlice1, resolvedSlice2] = ...
    happinessCheckBTD(fibersValues, colsValues, rowsValues, dom, pref, tech)
vsclFib = max(abs(fibersValues(:,1)));
vsclSli = max(abs(colsValues{1}(:)));
vsclSli = max(vsclSli, max(abs(rowsValues{1}(:))));
fiber1Data.hscale = norm(dom(5:6), inf);
fiber1Data.vscale = vsclFib;
fiber2Data.hscale = norm(dom(1:2), inf);
fiber2Data.vscale = vsclSli;
fiber3Data.hscale = norm(dom(3:4), inf);
fiber3Data.vscale = vsclSli;

fiber1Chebtech = tech.make(sum(fibersValues, 2), fiber1Data);
resolvedFibers  = happinessCheck(fiber1Chebtech, [], sum(fibersValues, 2), ...
    [], pref);

slice1Chebtech = tech.make(sum(colsValues{1}, 2), fiber2Data);
resolvedSlice1  = happinessCheck(slice1Chebtech, [], sum(colsValues{1}, 2),...
    [], pref);

slice2Chebtech = tech.make(sum(rowsValues{1}, 2), fiber3Data);
resolvedSlice2  = happinessCheck(slice2Chebtech, [], sum(rowsValues{1}, 2),...
    [], pref);

end

%%
function [op, dom, pref, vscaleBnd, fiberDim, vectorize, isEqui, ...
    fixedRank] = parseInputs(op, varargin)
vectorize = 0;
isEqui = 0;
isCoeffs = 0;
fixedRank = 0;
pref = chebfunpref();
fiberDim = [];
vscaleBnd = [];
dom = [-1 1 -1 1 -1 1];

% Preferences structure given?
isPref = find(cellfun(@(p) isa(p, 'chebfunpref'), varargin));
if ( any(isPref) )
    pref = varargin{isPref};
    varargin(isPref) = [];
end


if ( isa(op, 'char') )     % CHEBFUN3( CHAR )
    op = str2op(op);
end

for k = 1:length(varargin)
    if any(strcmpi(varargin{k}, {'trig', 'periodic'}))
        pref.tech = @trigtech;
    elseif strcmpi(varargin{k}, 'eps')
        pref.cheb3Prefs.chebfun3eps = varargin{k+1};
    elseif strcmpi(varargin{k}, 'rank') % rank is specified.
        fixedRank = varargin{k+1};
    elseif ( isnumeric(varargin{k}) )
        if ( numel(varargin{k}) == 6 ) % domain is specified.
            dom = varargin{k};   
        elseif ( numel(varargin{k}) == 3 ) % length is specified.
            if ( k > 1 && strcmpi(varargin{k-1}, 'rank') )
                % Rank is specified, not length. Don't confuse them.
                continue
            else
                % Interpret this as the user wants a fixed degree chebfun3 
                % on the domain DOM.
                len = varargin{k};
                [xx, yy, zz] = chebfun3.chebpts3(len(1), len(2), len(3), dom);
                op = op(xx, yy, zz);
            end
        end
    elseif strcmpi(varargin{k}, 'vscaleBnd')
        vscaleBnd = varargin{k+1};        
    elseif strcmpi(varargin{k}, {'fiberDim'})
        fiberDim = varargin{k+1};
    elseif any(strcmpi(varargin{k}, {'vectorize', 'vectorise'}))
        vectorize = true;
    elseif strcmpi(varargin{k}, 'coeffs')
        isCoeffs = 1;
    elseif strcmpi(varargin{k}, 'equi')
        isEqui = 1;
    end
end

if ( isCoeffs )
    op = chebfun3.coeffs2vals(op);
end

% If the vectorize flag is off, do we need to give user a warning?
if ( ~vectorize && ~isnumeric(op) ) % another check
    [vectorize, op] = vectorCheck(op, dom);
end

end

%% 
function [vectorize, op] = vectorCheck(op, dom)
% Check for cases like op = @(x,y,z) x*y^2*z

vectorize = false;
[xx, yy, zz] = ndgrid(dom(1:2), dom(3:4), dom(5:6));
try
    A = op(xx, yy, zz);
catch
    throwVectorWarning();
    vectorize = true;
    return
end

A = op(xx, yy, zz);
if ( any(isinf(A(:) ) ) )
    error('CHEBFUN:CHEBFUN3:constructor:inf', ...
        'Function returned INF when evaluated');
    elseif ( any(isnan(A(:)) ) )
        error('CHEBFUN:CHEBFUN3:constructor:nan', ...
            'Function returned NaN when evaluated');
end
if ( isscalar(A) )
    op = @(x,y,z) op(x,y,z) + 0*x + 0*y + 0*z;
end
end

%%
function throwVectorWarning()
warning('CHEBFUN:CHEBFUN3:constructor:vectorize',...
    ['Function did not correctly evaluate on an array.\n', ...
    'Turning on the ''vectorize'' flag. Did you intend this?\n', ...
    'Use the ''vectorize'' flag in the CHEBFUN3 constructor\n', ...
    'call to avoid this warning message.']);
end
%%
function op = str2op(op)
% OP = STR2OP(OP), finds independent variables in a string and returns an op
% handle than can be evaluated.

vars = symvar(op); % Independent variables
if ( numel(vars) > 3)
    error('CHEBFUN:CHEBFUN3:constructor:str2op:depvars', ...
        'Too many independent variables in string input.');
else
    op = eval(['@(' vars{1} ',' vars{2} ',' vars{3} ')' op]);
end

end
%%

function f = fixTheRank(f , fixedRank)
% Fix the rank of a CHEBFUN3. Used for calls to constructor with specified 
% rank.

if ( any(fixedRank) < 0 )
    error('CHEBFUN:CHEBFUN3:constructor:fixTheRank:negative', ...
        'Ranks should all be nonnegative.')
elseif ( all(fixedRank) )
    [r1, r2, r3] = rank(f);
    t1 = fixedRank(1);
    t2 = fixedRank(2);
    t3 = fixedRank(3);
    
    % What to do with cols?
    if ( r1 > t1 )
        % Truncate cols:
        f.cols = f.cols(:, 1:t1);
        f.core = f.core(1:t1, :, :);
        r1 = t1; % New size of core
    elseif ( r1 < t1 )
        % Pad cols with approprate number of zero cols:
        zCols = chebfun(0, f.cols.domain);
        for jj = r1 : t1 - 1
            f.cols = [f.cols zCols];
        end
        % Pad mode 1 of the core tensor with zeros:
        tempCore = zeros(t1, r2, r3);
        tempCore(1:r1, :, :) = f.core;
        f.core = tempCore;
        r1 = t1; % New size of core
    end
    
    % What to do with rows?
    if ( r2 > t2 )
        % Truncate rows:
        f.rows = f.rows(:,1:t2);
        f.core = f.core(:, 1:t2, :);
        r2 = t2; % New size of core
    elseif ( r2 < t2 )
        % Pad rows with approprate number of zero rows:
        zRows = chebfun(0, f.rows.domain);
        for jj = r2 : t2 - 1        
            f.rows = [f.rows zRows];
        end
        % Pad mode 2 of the core tensor with zeros:
        tempCore = zeros(r1, t2, r3);
        tempCore(:, 1:r2, :) = f.core;
        f.core = tempCore;
        r2 = t2; % New size of core
    end
    
    % What to do with tubes?
    if ( r3 > t3 )
        % Truncate tubes:
        f.tubes = f.tubes(:, 1:t3);
        f.core = f.core(:, :, 1:t3);
    elseif ( r3 < t3 )
        % Pad tubes with approprate number of zero tubes:
        zTubes = chebfun(0, f.tubes.domain);
        for jj = r3 : t3 - 1
            f.tubes = [f.tubes zTubes];
        end
        % Pad mode 3 of the core tensor with zeros:
        tempCore = zeros(r1, r2, t3);
        tempCore(:, :, 1:r3) = f.core;
        f.core = tempCore;
    end
end

end

%%
function f = constructFromDouble(op, dom, pref, isEqui)

pseudoLevel = pref.cheb3Prefs.chebfun3eps;

f = chebfun3();
if ( ~isEqui && numel(op) == 1 )
    f = constructor(f, @(x,y,z) op + 0*x, dom);
    return;
end

if numel(size(op)) == 3 % We have a tensor of values.
    % Calculate a tolerance and find numerical rank to this tolerance:
    % The tolerance assumes the samples are generated by NDGRID from a 
    % function. It depends on the size of the sample tensor, hscale of 
    % domain, vscale of the samples, condition number of the function, and 
    % the accuracy target in chebfun3 preferences.
    % N.B. We cannot detect if MESHGRID was used to generate values unless
    % we know also the (x,y,z) points used to generate those values.
    % If we knew beforehand that ALL users WILL generate their tensor of 
    % values ONLY from meshgrid pts, all we need is to say 
    % vals = permute(vals,[2 1 3]); to generate a tensor corresponding to 
    % ''meshgrid'', in which case a copy of ``points3D`` should also be 
    % used accordingly. op = permute(op,[2 1 3]);
    
    if ( ~isEqui )
        [xx, yy, zz] = points3D(size(op,1), size(op,2), size(op,3), dom, pref);
    else
        % Equispaced points from ndgrid, not meshgrid!
        x = linspace(dom(1), dom(2), size(op, 1));
        y = linspace(dom(3), dom(4), size(op, 2));
        z = linspace(dom(5), dom(6), size(op, 3));
        [xx, yy, zz] = ndgrid(x, y, z);
    end
    [relTol, absTol] = getTol3D(xx, yy, zz, op, max(size(op)), dom, ...
        pseudoLevel, []);
    pref.chebfuneps = relTol;

    
    % Perform 3D ACA with complete pivoting:
    fiberDim = 3;
    factor = 0;
    [colsValues, rowsValues, pivotVals2D, ignore, tubesValues, ...
        pivotVals3D, ignore, ignore, ignore] = completeACA3D(op, ...
        fiberDim, absTol, factor, dom, pref);
    
    sepRank = numel(pivotVals3D); % first separation rank
    for k=1:sepRank
        diagValues2D{k} = diag(1./pivotVals2D{k});
    end
    
    % BTD ---> Tucker compression:
    dom2D = dom(1:4); % The first 4 entries correspond to slices
    [core, colsValues, rowsValues] = btd2tucker(colsValues, rowsValues, ...
        diagValues2D, pivotVals3D, dom2D, pref, absTol);
    
    % Construct a CHEBFUN3 object and call a sampleTest.
    if ( ~isEqui )
        f.cols = chebfun(colsValues, dom(1:2), pref);
        f.rows = chebfun(rowsValues, dom(3:4), pref);
        f.tubes = chebfun(tubesValues, dom(5:6), pref);
    else
        f.cols = chebfun(colsValues, dom(1:2), 'equi', pref);
        f.rows = chebfun(rowsValues, dom(3:4), 'equi', pref);
        f.tubes = chebfun(tubesValues, dom(5:6), 'equi', pref);
    end
    
    % TODO: Do we want to simplify when constructing from doubles?
    % The following causes the 2nd test in test_construnctorsyntax() fail.
%     f.cols = simplify(f.cols, pref.chebfuneps, 'globaltol');
%     f.rows = simplify(f.rows, pref.chebfuneps, 'globaltol');
%     f.tubes = simplify(f.tubes, pref.chebfuneps, 'globaltol');
    
    f.core = core;
    f.domain = dom;
    return;
end

end

%%
function [colsBtd, rowsBtd, pivotValues2D, pivotIndices2D, fibers, ...
    pivotValues3D, pivotIndices3D, ifail3D, ifail2D] = completeACA3D(A, ...
    fiberDim, tol, factor, dom, pref)
%   Non-adaptive (fixed-size) MACA, i.e., a 3D analogue of Gaussian 
%   elimination with complete pivoting.
%
%   INPUTS:     A:        A given tensor of function values at 3D chebpts.
%
%               fiberDim: Dimension to be used for 1st separation of fibers
%                         and slices.
%
%               tol:      A given tolerance on the magnitude of the pivots.
%
%               factor:   The ratio between the width of A and the number 
%                         of iterations (rank) allowed in Gaussian elimination.
%
%  OUTPUTS:     colsBtd: A cell-array containing skeleton columns of slices 
%                        in block term decomposition.
%
%               rowsBtd: A cell-array containing skeleton rows of slices in
%                        block term decomposition.
%
%               pivotValues2D: A cell array containing the values of pivots
%                              in 2D ACAs.
%
%               pivotIndices2D: A cell array containing indices i.e.,
%                               locations of of 2D pivot points.
%
%                       fibers: A matrix of size n1 x ITER. Each of its columns
%                             contains the values of the updated (residue) 
%                             tensor at the pivot fiber.
%
%               pivotValues3D: A row vector containing the values of the 
%                             pivot entries during the iterations. 
%
%               pivotIndices3D: A matrix of size rank x 3, where rank = iter. 
%                         Each of its rows contain the index of one 3D pivotValues.
%
%                      ifail: We fail if iter >= (width/factor).

% Developer Note: The output of this code should satisfy the following 
% slice decomposition:
%       AA \approx temp, 
% where AA is a copy of A from input, and temp is computed as follows:
%   temp = zeros(size(A2)); 
%   for i = 1:3, 
%    temp = temp + chebfun3.outerProd(slices(:,:,i),fibers(:,i)./pivotValues3D(i));
%   end
%   norm(AA(:) - temp(:))
%
% An analogous BTD decomposotion should also hold.

pseudoLevel = pref.cheb3Prefs.chebfun3eps;

% Set up output variables.
[n1, n2, n3] = size(A);
if fiberDim == 1
    width = min(n1, n2*n3);    % Use to tell us how many pivots we can take
elseif fiberDim == 2
    width = min(n2, n1*n3);   
else
   width = min(n3, n1*n2);
end
pivotValues3D = zeros(1);      % Store an unknown number of Pivot values
pivotIndices3D = zeros(1, 3);  % Store (col, row, tube) = entries of pivot location
ifail3D = 1;                   % Assume we fail in 3D ACA
ifail2D = 1;                   % Assume we also fail in the 2D ACAs
globalTol = [];
sliceDim = [1 2];              % See Developer note in the following.

% Main algorithm
iter = 0;                  % Count number of interpolated rows/slices.
[infNorm, ind] = max(abs(reshape(A, numel(A), 1))); % Complete pivoting
[col, row, tube] = ind2sub(size(A), ind);

scl = infNorm;
% If the function is the zero function.
if ( scl == 0 )
    pivotValues3D = 0;
    fibers = 0;
    colsBtd{1} = 0;
    rowsBtd{1} = 0;
    pivotValues2D = 0;
    pivotIndices2D = [0 0];
    ifail3D = 0;
    ifail2D = 0;
else
    fibers(:,1) = zeros(size(A, 3),1);
    colsBtd{1} = zeros(size(A, sliceDim(1)),1);
    rowsBtd{1} = zeros(size(A, sliceDim(2)),1);
    pivotValues2D{1} = 0;
    pivotIndices2D{1} = [0 0];
    slices(:,:,1) = zeros(size(A,sliceDim(1)), size(A, sliceDim(2)), 1);
end
dom2D = dom(1:4);

while ( ( infNorm > tol ) && ( iter < width / factor) ...
        && ( iter < width ) )
    fibers(:, iter+1) = A(col, row, :);  % Extract skeleton tubes. Each
    % column in "fibers" is N3 x 1.
    
    slices(:,:,iter+1) = A(:,:,tube);    % Extract skeleton slices to be 
    % decomposed further. Each slice in "slices" is N1 x N2.
    
    % Developer Note: As the above lines show, we are always separating the
    % last variable from the first two. The point is that the function 
    % handle has already been permuted outside this subroutine and 
    % therefore the tensor A here contains values of the permuted function.
    % In this sense, here we are in essense separating the variable chosen 
    % by the dimension clustering step, and not necessarily the last 
    % variable z.
    
    PivVal3D = A(col, row, tube);        % = f(X(col), Y(row), Z(tube)) in
                                         % the 1st iteration.
    
    % Use the first slice to compute globalTol for 2D ACAs applied to all 
    % slices.
    if iter == 0
        [xx2D, yy2D] = points2D(n1, n2, dom2D, pref); % ndgrid is used.
        globalTol = GetTol2Dv2(xx2D, yy2D, slices(:,:,1), dom2D, pseudoLevel);
    end
    
    % Apply 2D ACA to each slice to form columns and rows in block term
    % decomposition:
    [colsBtd{iter+1}, pivotValues2D{iter+1}, rowsBtd{iter+1}, ...
        pivotIndices2D{iter+1}, ifail2DIter] = ...
        chebfun2ACAv2(slices(:, :, iter+1), globalTol, factor);
    
    % Developer Note: Since we use globalTol for slices after 1st 
    % iteration, it might be that these 2D ACA's don't fail, while with a 
    % localTol they would fail. So, it is auaully the 1st slice which shows
    % whether or not we got the 2D rank right. So, just use that one:
    if iter == 0
        ifail2D = ifail2DIter;
    end
    
    % Update the tensor, i.e., compute the residual tensor:
    A = A - chebfun3.outerProd(colsBtd{iter+1} * ...
        (diag(1./(pivotValues2D{iter+1}))) * (rowsBtd{iter+1}.'), ...
        fibers(:,iter+1) ./ PivVal3D);
    % Equivalently, we have: 
    % A = A - chebfun3.outerProd(slices(:,:,iter+1),fibers(:,iter+1)./PivVal3D);
    
    % Keep track of progress in 3D ACA:
    iter = iter + 1;              % One more fiber and slice are removed from A
    pivotValues3D(iter) = PivVal3D;           % Store pivot value in 3D ACA
    pivotIndices3D(iter, :)=[col row tube];    % Store pivot location in 3D ACA
    
    % Find next 3D pivot value and its location:
    [infNorm, ind] = max(abs(A(:)));
    [col, row, tube] = ind2sub(size(A), ind);
end

if ( iter>0 && all(pivotValues2D{iter} == 0) )    
    % If the last 2D pivot was zero, remove it
    colsBtd=colsBtd(1:iter-1);
    rowsBtd=rowsBtd(1:iter-1);
    pivotValues2D = pivotValues2D(1:iter-1);
    pivotIndices2D = pivotIndices2D(1:iter-1);
    fibers = fibers(:,1:iter-1);
    pivotValues3D = pivotValues3D(1:iter-1);
    pivotIndices3D=pivotIndices3D(1:iter-1,:);
    infNorm = 0; % If the last 2D pivot was zero, infNorm will be NaN and 
    % it makes the next statement result in ifail3D \neq 0; So, put it 0 to
    % get ifail3D = 0;
end

if ( infNorm <= tol )
    ifail3D = 0;                               % We didn't fail in 3D ACA.
end
if ( iter > width/factor )
    ifail3D = 1;                               % We did fail in 3D ACA.
end

end
%%
