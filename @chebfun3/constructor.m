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
[op, dom, pref, fiberDim, vectorize, isEqui, fixedRank] = ...
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

%% Dimension clustering
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
    out = tech.tensorGrid([grid, grid, grid], dom);
    xx = out{1};
    yy = out{2};
    zz = out{3};
    grid2D = grid;
    vals = evaluate(op, xx, yy, zz, vectorize);
    % We have vals(i,j,k) = op(X(i), Y(j), Z(k)), where                 (*)
    %                                           X = chebpts(m,[a,b]);
    %                                           Y = chebpts(n,[c,d]);
    %                                           Z = chebpts(p,[e,g]);
    %                                 [xx, yy, zz] = ndgrid(X, Y, Z);
    % and ndgrid is used here in ``tensorGrid``.
    % If we have used meshgrid, instead of ndgrid, in ``tensorGrid``, then 
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
    [relTol, absTol] = getTol3D(xx, yy, zz, vals, grid, dom, pseudoLevel);
    pref.chebfuneps = relTol; % tolerance to be used in happinessCheck.
    
    %% PHASE 1: %%
    % Apply 3D ACA to tensor of values
    [colsValues, rowsValues, pivotVals2D, pivotIndices2D, fibersValues, ...
        pivotVals3D, pivotIndices3D, iFail3D, iFail2D] = completeACA3D(...
        vals, absTol, factor, dom, pref);
    
    while ( (iFail3D || iFail2D) && grid < maxSamplePhase1 )
        % Refine sampling on tensor grid:
        if ( iFail3D )
            grid = gridRefinePhase1(grid, pref);
            out = tech.tensorGrid([grid, grid, grid], dom);
            xx = out{1};
            yy = out{2};
            zz = out{3};
        elseif ( iFail2D )
            grid2D = gridRefinePhase1(grid2D, pref);
            % TODO: Do we need a bound on grid2D?
            out = tech.tensorGrid([grid2D, grid2D, grid], dom);
            xx = out{1};
            yy = out{2};
            zz = out{3};
        end
        
        vals = evaluate(op, xx, yy, zz, vectorize); % Resample
        % New tolerance:
        [relTol, absTol] = getTol3D(xx, yy, zz, vals, grid, dom, ...
            pseudoLevel);
        pref.chebfuneps = relTol;
        
        % New 3D ACA:
        [colsValues, rowsValues, pivotVals2D, pivotIndices2D, ...
            fibersValues, pivotVals3D, pivotIndices3D, iFail3D, iFail2D] ...
            = completeACA3D(vals, absTol, factor, dom, pref);
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
        diagValues2D = cell(sepRank, 1);
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
            Z = tech.tensorGrid(p, dom(5:6));
            fibersValues = zeros(p, sepRank);
            for k=1:sepRank
                [xx, yy, zz] = ndgrid(pivPos3D(k, 1), pivPos3D(k, 2), Z);
                fibersValues(:,k) = squeeze(evaluate(op, xx, yy, zz, vectorize));
            end
            
            % Find location of pivots on new grid  (using nesting property).
            PI3D(:, 3) = nesting(PI3D(:, 3)); % The first column of PP contains 
            % z indices and should therefore be updated.
        else
            fibersValues = zeros(p, sepRank);
            Z = tech.tensorGrid(p, dom(5:6));
            for k=1:sepRank
                [xx, yy, zz] = ndgrid(pivPos3D(k, 1), pivPos3D(k, 2), Z);
                fibersValues(:, k) = squeeze(evaluate(op, xx, yy, zz, vectorize));
            end
        end
        
        if ( ~resolvedCols && resolvedRows )
            % Double sampling along the column slices only, i.e., in x. 
            [m, nesting] = gridRefinePhase2(m, pref); 
            X = tech.tensorGrid(m, dom(1:2));
            Y = tech.tensorGrid(n, dom(3:4));
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
            X = tech.tensorGrid(m, dom(1:2));
            Y = tech.tensorGrid(n, dom(3:4));
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
            X = tech.tensorGrid(m, dom(1:2));
            Y = tech.tensorGrid(n, dom(3:4));
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
        colsData.hscale = norm(dom(1:2), inf);
        colsData.vscale = max(abs(colsValues{1}(:)));        
        if ( ~resolvedCols )
            colsChebtech = tech.make(colsValues{1}, colsData);
            colsChebtech.coeffs = sum(abs(colsChebtech.coeffs), 2);
            resolvedCols  = happinessCheck(colsChebtech, [], ...
                colsChebtech.coeffs, [], pref);
        end
        
        rowsData.hscale = norm(dom(3:4), inf);
        rowsData.vscale = max(abs(rowsValues{1}(:)));
        if ( ~resolvedRows )
            rowsChebtech = tech.make(rowsValues{1}, rowsData);
            rowsChebtech.coeffs = sum(abs(rowsChebtech.coeffs), 2);
            resolvedRows  = happinessCheck(rowsChebtech, [], ...
                rowsChebtech.coeffs, [], pref);
        end

        tubesData.hscale = norm(dom(5:6), inf);
        tubesData.vscale = max(abs(fibersValues(:)));
        if ( ~resolvedTubs )          
            fiber3Chebtech = tech.make(fibersValues, tubesData);
            fiber3Chebtech.coeffs = sum(abs(fiber3Chebtech.coeffs), 2);
            resolvedTubs  = happinessCheck(fiber3Chebtech, [], ...
                fiber3Chebtech.coeffs, [], pref);
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
        [core, colsValues, rowsValues] = btd2tucker(colsValues, ...
            rowsValues, diagValues2D, pivotVals3D, absTol);
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
    rowsValues, diagValues2D, pivotVals3D, absTol)
allCols = []; 
allRows = []; 
allDiags = [];
nn = numel(pivotVals3D);
sizeIndex = zeros(nn, 1);

for kkk = 1:nn
    sizeIndex(kkk) = size(rowsValues{kkk}, 2);
    allCols = [allCols, colsValues{kkk}];
    allRows = [allRows, rowsValues{kkk}];
    allDiags = [allDiags; diag(diagValues2D{kkk})];
end
sizeIndex = cumsum([0; sizeIndex]);
        
% Compress allCols and allRows:
if ( size(allCols, 2) > 1 )
    [Su, ~, Vu, colsTucker] = completeACA2D(allCols, absTol, 0);
    % factor = 0, because we want the ACA to be applied even if op is not 
    % low-rank. In contrast to Chebfun2, we now have 
    % allCols = colsTucker * diag(1./Su) * Vu'.
    % Moreover, we use the same absTol to chop columns, the same tolerance
    % used to chop fibers so that ranks are less inconsistent for symmetric
    % functions.
    Vu = Vu*diag(1./Su);
else
    colsTucker = allCols; 
    Vu = 1;
end

if ( size(allRows, 2) > 1 )
    [Sv, ~, Vv, rowsTucker] = completeACA2D(allRows, absTol, 0);
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
function absTol = GetTol2D(xx, yy, vals, dom, pseudoLevel)
% GETTOL2D   Calculate a tolerance for the Chebfun2 constructor.
%
%  This is the 2D analogue of the tolerance employed in the chebtech
%  constructors. It is based on a finite difference approximation to the
%  gradient, the size of the approximation domain, the internal working
%  tolerance, and an arbitrary (4/5) exponent. 

[m, n] = size(vals); 
grid = max(m, n);
relTol = grid^(4/5) * pseudoLevel; % this should be vscale and hscale invariant

% Remove some edge values so that df_dx and df_dy have the same size. 
% xx is generated by ndgrid, i.e., xx changes in the first mode:
dfdx = diff(vals(:, 1:n-1), 1, 1) ./ diff(xx(:, 1:n-1), 1, 1);
% yy is generated by ndgrid, i.e., yy changes row-wise (2nd mode):
dfdy = diff(vals(1:m-1, :), 1, 2) ./ diff(yy(1:m-1, :), 1, 2);
gradNorms = [max(abs(dfdx(:))), max(abs(dfdy(:)))];
% A vector of gradient information over the domain.
if ( isempty(gradNorms) )
    % This happens if the input in not a trivariate function.
    gradNorms = 1;
end

vscale = max(abs(vals(:)));
domDiff = [diff(dom(1:2)) diff(dom(3:4))];
absTol = max(max(gradNorms.*domDiff), vscale) * relTol;
end

%%
function [col, pivotVals, row, pivotLoc, ifail2D] = chebfun2ACA(op, ...
    tol, factor)
% Perform GE with complete pivoting:

if ( factor ~= 0 )
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
    pseudoLevel)

relTol = 2*grid^(4/5) * pseudoLevel; % this should be vscale and hscale invariant
vscale = max(abs(vals(:)));
[m,n,p] = size(vals);
% Remove some edge values so that df_dx, df_dy and df_dz have the same size. 
% xx changes in the first mode:
df_dx = diff(vals(:, 1:n-1, 1:p-1), 1, 1) ./ diff(xx(:, 1:n-1, 1:p-1), 1, 1);
% yy changes row-wise (2nd mode):
df_dy = diff(vals(1:m-1, :, 1:p-1), 1, 2) ./ diff(yy(1:m-1, :, 1:p-1), 1, 2);
% zz changes tube-wise (3rd mode):
df_dz = diff(vals(1:m-1, 1:n-1, :), 1, 3) ./ diff(zz(1:m-1, 1:n-1, :), 1, 3);
gradNorms = [max(abs(df_dx(:))), max(abs(df_dy(:))), max(abs(df_dz(:)))];
% A vector of gradient information over the domain.
if ( isempty(gradNorms) )
    % This happens if the input in not a trivariate function in which case
    % we basically disable using gradient information:
    gradNorms = 1;
end
domDiff = [diff(dom(1:2)) diff(dom(3:4)) diff(dom(5:6))];
absTol = max(max(gradNorms.*domDiff), vscale) * relTol;
% absTol should depend on the vscale of the function while it also uses
% derivative information to prevent issues like the one mentioned in
% https://github.com/chebfun/chebfun/issues/1491.
    
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
tech = pref.tech();
out = tech.tensorGrid([grid, grid, grid], dom);
xx = out{1};
yy = out{2};
zz = out{3};
vals = evaluate(op, xx, yy, zz, vectorize);

% Using SVD and 1D Chebfun:
F1 = chebfun3.unfold(vals, 1);
F2 = chebfun3.unfold(vals, 2);
F3 = chebfun3.unfold(vals, 3);
rX = rank(F1.'); % Transpose to have longer columns
rY = rank(F2.'); % and therefore being faster
rZ = rank(F3.'); % in MATLAB.
r = [rX, rY, rZ];
[~, ind] = find(r==min(r));
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
    [~, index] = max([len(ind(1)), len(ind(2))]);
    fiberDim = ind(index);
else % numel(ind) = 3
    [~, fiberDim] = max(len);
end

end

%%
function [resolvedFibers, resolvedSlice1, resolvedSlice2] = ...
    happinessCheckBTD(fibersValues, colsValues, rowsValues, dom, pref, tech)
% Happiness-check for the BTD: It uses all the fibers and columns and rows 
% in the first term of the BTD.
vsclFib = max(abs(fibersValues(:,1)));
vsclSli = max(abs(colsValues{1}(:)));
vsclSli = max(vsclSli, max(abs(rowsValues{1}(:))));
fiber1Data.hscale = norm(dom(5:6), inf);
fiber1Data.vscale = vsclFib;
fiber2Data.hscale = norm(dom(1:2), inf);
fiber2Data.vscale = vsclSli;
fiber3Data.hscale = norm(dom(3:4), inf);
fiber3Data.vscale = vsclSli;

% Convert to coefficients:
fiber1Chebtech = tech.make(fibersValues, fiber1Data); 
% Add absolute value of all the coefficients:
fiber1Chebtech.coeffs = sum(abs(fiber1Chebtech.coeffs), 2);
% Check for happiness using the corresponding tech:
resolvedFibers  = happinessCheck(fiber1Chebtech, [], ...
    fiber1Chebtech.coeffs, [], pref);

slice1Chebtech = tech.make(colsValues{1}, fiber2Data);
slice1Chebtech.coeffs = sum(abs(slice1Chebtech.coeffs), 2);
resolvedSlice1  = happinessCheck(slice1Chebtech, [], ...
    slice1Chebtech.coeffs, [], pref);
    
slice2Chebtech = tech.make(rowsValues{1}, fiber3Data);
slice2Chebtech.coeffs = sum(abs(slice2Chebtech.coeffs), 2);
resolvedSlice2  = happinessCheck(slice2Chebtech, [], ...
    slice2Chebtech.coeffs, [], pref);

end

%%
function [op, dom, pref, fiberDim, vectorize, isEqui, ...
    fixedRank] = parseInputs(op, varargin)
vectorize = 0;
isEqui = 0;
isCoeffs = 0;
fixedRank = 0;
pref = chebfunpref();
fiberDim = [];
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
                [xx, yy, zz] = chebpts3(len(1), len(2), len(3), dom);
                op = op(xx, yy, zz);
            end
        end
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
    A = feval(op, xx, yy, zz);
catch
    throwVectorWarning();
    vectorize = true;
    return
end

A = feval(op, xx, yy, zz);
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
% OP = STR2OP(OP), finds independent variables in a string and returns an 
% op handle than can be evaluated.

vars = symvar(op);        % Independent variables
numVars = numel(vars);
if ( numVars == 0 )
    op = @(x,y,z) eval (op);
    
elseif ( numVars == 1 )
    op = eval(['@(' vars{1} ', myVarBeta, myVarGamma)' op]);
    
elseif ( numVars == 2 )
    op = eval(['@(' vars{1} ',' vars{2} ', myVarGamma)' op]);

elseif ( numVars == 3 )
    op = eval(['@(' vars{1} ',' vars{2} ',' vars{3} ')' op]);

else
    error('CHEBFUN:CHEBFUN3:constructor:str2op:depvars', ...
        'Too many independent variables in string input.');
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
% The input is a discrete tensor of values.
pseudoLevel = pref.cheb3Prefs.chebfun3eps;
tech = pref.tech();

f = chebfun3();
if ( ~isEqui && numel(op) == 1 )
    f = constructor(f, @(x,y,z) op + 0*x, dom);
    return;
end

% N.B. We cannot detect if MESHGRID was used to generate values unless we 
% know also the (x,y,z) points used to generate those values. If we knew 
% beforehand that ALL users WILL generate their tensor of values ONLY from 
% meshgrid pts, all we need is to say vals = permute(vals,[2 1 3]); to 
% generate a tensor corresponding to ''meshgrid'', in which case a copy of 
% ``tensorGrid`` should also be used accordingly. op = permute(op,[2 1 3]);
if ( ~isEqui )
    m = size(op, 1);
    n = size(op, 2);
    p = size(op, 3);
    out = tech.tensorGrid([m, n, p], dom);
    xx = out{1};
    yy = out{2};
    zz = out{3};
else
    % Equispaced points from ndgrid, not meshgrid!
    x = linspace(dom(1), dom(2), size(op, 1));
    y = linspace(dom(3), dom(4), size(op, 2));
    z = linspace(dom(5), dom(6), size(op, 3));
    [xx, yy, zz] = ndgrid(x, y, z);
end

% Calculate a tolerance and find numerical rank to this tolerance: The 
% tolerance assumes the samples are generated by NDGRID from a function. 
% It depends on the size of the sample tensor, hscale of domain, vscale of 
% the samples, condition number of the function, and the accuracy target in
% chebfun3 preferences.
[relTol, absTol] = getTol3D(xx, yy, zz, op, max(size(op)), dom, pseudoLevel);
pref.chebfuneps = relTol;

% Perform 3D ACA with complete pivoting:
factor = 0;
[colsValues, rowsValues, pivotVals2D, ~, tubesValues, pivotVals3D, ~, ~, ...
    ~] = completeACA3D(op, absTol, factor, dom, pref);

sepRank = numel(pivotVals3D); % first separation rank
diagValues2D = cell(sepRank, 1);
for k=1:sepRank
        diagValues2D{k} = diag(1./pivotVals2D{k});
end

% BTD ---> Tucker compression:
[core, colsValues, rowsValues] = btd2tucker(colsValues, rowsValues, ...
    diagValues2D, pivotVals3D, absTol);

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
f.core = core;
f.domain = dom;
return

end

%%
function [colsBtd, rowsBtd, pivotValues2D, pivotIndices2D, fibers, ...
    pivotValues3D, pivotIndices3D, ifail3D, ifail2D] = completeACA3D(A, ...
    tol, factor, dom, pref)
%   Non-adaptive (fixed-size) MACA, i.e., a 3D analogue of Gaussian 
%   elimination with complete pivoting.
%
%   INPUTS:     A:        A given tensor of function values at 3D chebpts.
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
tech = pref.tech();

% Set up output variables.
[n1, n2, n3] = size(A);
width = min(n3, n1*n2);        % Use to tell us how many pivots we can take
                               % See Developer note in the following.
pivotValues3D = zeros(1);      % Store an unknown number of Pivot values
pivotIndices3D = zeros(1, 3);  % Store (col, row, tube) = entries of pivot location
ifail3D = 1;                   % Assume we fail in 3D ACA
ifail2D = 1;                   % Assume we also fail in the 2D ACAs
globalTol = [];
sliceDim = [1 2];              % See Developer note in the following.

% Main algorithm
iter = 0;                      % Count number of interpolated rows/slices.
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
    fibers(:,1) = zeros(size(A, 3), 1);
    colsBtd{1} = zeros(size(A, sliceDim(1)), 1);
    rowsBtd{1} = zeros(size(A, sliceDim(2)), 1);
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
    if ( iter == 0 )
        out = tech.tensorGrid([n1, n2], dom2D);
        xx2D = out{1};
        yy2D = out{2};
        globalTol = GetTol2D(xx2D, yy2D, slices(:,:,1), dom2D, pseudoLevel);
    end
    
    % Apply 2D ACA to each slice to form columns and rows in block term
    % decomposition:
    [colsBtd{iter+1}, pivotValues2D{iter+1}, rowsBtd{iter+1}, ...
        pivotIndices2D{iter+1}, ifail2DIter] = ...
        chebfun2ACA(slices(:, :, iter+1), globalTol, factor);
    
    % Developer Note: Since we use globalTol for slices after 1st 
    % iteration, it might be that these 2D ACA's don't fail, while with a 
    % localTol they would fail. So, it is auaully the 1st slice which shows
    % whether or not we got the 2D rank right. So, just use that one:
    if ( iter == 0 )
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

if ( ( iter > 0 ) && ( all(pivotValues2D{iter} == 0) ) )
    % If the last 2D pivot was zero, remove it
    colsBtd=colsBtd(1:iter-1);
    rowsBtd=rowsBtd(1:iter-1);
    pivotValues2D = pivotValues2D(1:iter-1);
    pivotIndices2D = pivotIndices2D(1:iter-1);
    fibers = fibers(:, 1:iter-1);
    pivotValues3D = pivotValues3D(1:iter-1);
    pivotIndices3D=pivotIndices3D(1:iter-1,:);
    infNorm = 0; % If the last 2D pivot was zero, infNorm will be NaN and 
    % it makes the next statement result in ifail3D \neq 0; So, put it 0 to
    % get ifail3D = 0;
end

if ( infNorm <= tol )
    ifail3D = 0;                               % We didn't fail in 3D ACA.
    if ( iter == 0 )
        ifail2D = 0;    
    end
end
if ( iter > width/factor )
    ifail3D = 1;                               % We did fail in 3D ACA.
end

end
%%