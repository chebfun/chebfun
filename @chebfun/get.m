function out = get(f, prop)
%GET   GET method for the CHEBFUN class
%   P = GET(F, PROP) returns the property P specified in the string PROP from
%   the CHEBFUN F. Valid entries for the string PROP are:
%       'domain'         - The domain of definition of F.
%       'ends'
%       'funs'           - The piecewise smooth components of F.
%       'vscale'         - Vertical scale of F.
%       'vscale-local'   - Local vertical scales of F.
%       'hscale'         - Horizontal scale of F.
%       'hscale-local'   - Local horizontal scales of F.
%       'ishappy'        - Is F happy?
%       'epslevel'       - Approximate accuracy estimate of F.
%       'epslevel-local' - Approximate accuracy estimate of F's components.
%       'lval'           - Value(s) of F at left-hand side of domain.
%       'rval'           - Value(s) of F at right-hand side of domain.
%       'lval-local      - Value(s) of F's FUNs at left sides of their domains.
%       'rval-local'     - Value(s) of F's FUNs at right sides of their domains.
%       'exps'           - Exponents in a CHEBFUN, a two vector.
%       'exponents'
%       'deltas'         - Return the locations and magnitude of delta functions
%                          as a single matrix with the first row corresponding
%                          to the locations of the delta functions, and the
%                          magnitudes appended below the first row.
%       'imps'           - Same as DELTAS, supported for backward compatibility.
%   The following are also supported for backward compatibility, and really only
%   make sense when the CHEBFUN is represented by a CHEBTECH-type object. Note
%   that in these cases a cell is always returned, even if the Chebfun has only
%   a single FUN.
%       'points'         - The Chebyshev grid used to represent F.
%       'values'         - The values of the CHEBFUN on the grid above.
%       'coeffs'         - The corresponding Chebyshev coefficients.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

switch prop
    case {'domain', 'ends'}
        out = domain(f);
    case 'vscale'
        out = vscale(f);
    case 'hscale'
        out = hscale(f);
    case 'epslevel'
        out = epslevel(f);
    case 'ishappy'
        out = ishappy(f);
    case fieldnames(f)
        out = cell(1, numel(f));
        for k = 1:1:numel(f)
            out{k} = f(k).(prop);
        end

        if ( numel(out) == 1 )
            out = out{1};
        end
    case 'lval'
        dom = domain(f);
        out = feval(f, dom(1));
    case 'rval'
        dom = domain(f);
        out = feval(f, dom(end));
    case 'lval-local'
        out = getSimpleNumericLocalProp(f, 'lval', true);
    case 'rval-local'
        out = getSimpleNumericLocalProp(f, 'rval', true);
    case 'vscale-local'
        out = getSimpleNumericLocalProp(f, 'vscale', true);
    case 'hscale-local'
        out = getSimpleNumericLocalProp(f, 'hscale', true);
    case 'epslevel-local'
        out = getSimpleNumericLocalProp(f, 'epslevel', true);
    case {'values', 'coeffs', 'points'}
        out = getSimpleNumericLocalProp(f, prop, true);
    case {'exps', 'exponents'}
        out = getArbitraryLocalProp(f, 'exponents', true);
        if ( isvector(out) )
            out = cell2mat(out);
        end
    case 'deltas'
        out = getArbitraryLocalProp(f, 'deltas', false);
        for k = 1:numColumns(f)
            colDeltaData = out{k}.';
            colDeltaData = horzcat(colDeltaData{:});
            if ( ~isempty(colDeltaData) )
                [mag, loc] = deltafun.mergeColumns(colDeltaData(2:end, :), ...
                    colDeltaData(1, :));
                out{k} = [loc ; mag];
            else
                out{k} = [];
            end
        end

        if ( numel(out) == 1 )
            out = out{1};
        end
    case 'deltas-local'
        out = getArbitraryLocalProp(f, 'deltas', true);
    case 'imps'
        warning('CHEBFUN:get:imps', ...
            '''imps'' property is deprecated.  Use ''deltas'' instead.');
        out = get(f, 'deltas');
    otherwise
        error('CHEBFUN:get:badProp', ...
            '''%s'' is not a recognized CHEBFUN property name.', prop);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = getSimpleNumericLocalProp(f, prop, doSimplify)
%GETSIMPLENUMERICLOCALPROP   Get simple numeric FUN properties.
%   OUT = GETSIMPLENUMERICLOCALPROP(F, PROP, DOSIMPLIFY) retrieves the values
%   of the "simple numeric local" property PROP from each of the FUNs of F and
%   assembles them in some sensible way.  A property is a "simple numeric
%   local" property if (1) it is numeric, (2) it is a FUN property, and (3)
%   either one of the following conditions holds:
%
%     A.  If F has N columns, then the property also has N columns, and the
%     values in column K of the property are associated with column K of F.
%
%     B.  The property consists of a single column, the values of which are
%     associated with _all_ columns of F.
%
%   If DOSIMPLIFY is TRUE, the function will try to represent the property in a
%   simpler form.  The hierarchy of "simplicity", from least simplest to
%   simplest, is as follows:
%
%     1.  Cell array of cell arrays.  out{j}{k} gives the value of the property
%     for FUN j of column k.
%
%     2.  2D cell array.  The output can be simplified to this if all columns
%     of F have the same number of FUNs.  out{j, k} gives the value of the
%     property for FUN j of column k.
%
%     3.  Numeric matrix.  If the output can be simplified to a 2D cell array,
%     we can sometimes convert it further into a numeric matrix if the result
%     can be unambiguously interpreted.  See CANCONVERT2DCELLARRAYTONUMERIC.
%
%   The reason this funciton exists is that simple numeric local properties can
%   be extracted from array-valued CHEBFUNs without needing to convert them to
%   quasimatrices first, which can potentially save a lot in computational
%   overhead costs.

% Get properties for each fun and assemble into cell arrays.
if (numel(f) == 1)  % Single column or array-valued CHEBFUN.
    out = getSimpleNumericLocalPropArrayValued(f, prop);
else                % Quasimatrix / array of CHEBFUNs.
    out = getArbitraryLocalPropQuasimatrix(f, prop);
end

% Try to convert from a cell array of cell arrays to a 2D cell array.
if ( doSimplify && canConvertCellOfCellsTo2DCellArray(out) )
    out = horzcat(out{:});

    % Try to go one step further and convert to a numeric matrix.
    if ( canConvert2DCellArrayToNumeric(out) )
        out = cell2mat(out);
    end
end

if ( f(1).isTransposed )
    out = out.';
end

end

function out = getArbitraryLocalProp(f, prop, doSimplify)
%GETSIMPLENUMERICLOCALPROP   Get simple numeric FUN properties.
%   OUT = GETSIMPLENUMERICLOCALPROP(F, PROP, DOSIMPLIFY) retrieves the values
%   of the property PROP from each of the FUNs of F and %   assembles them in
%   some sensible way.
%
%   If DOSIMPLIFY is TRUE, the function will try to represent the property in a
%   simpler form.  The hierarchy of "simplicity", from least simplest to
%   simplest, is as follows:
%
%     1.  Cell array of cell arrays.  out{j}{k} gives the value of the property
%     for FUN j of column k.
%
%     2.  2D cell array.  The output can be simplified to this if all columns
%     of F have the same number of FUNs.  out{j, k} gives the value of the
%     property for FUN j of column k.
%
%     3.  Single cell contents.  If the output can be simplified to a 2D cell
%     array, with only one cell, we just return the contents of that cell.

% The only way to correctly support access to general properties of
% array-valued CHEBFUNs is to convert them to quasimatrices first.  (Example:
% If f is an array-valued CHEBFUN with n columns and has a tech with a property
% that is a 1 x n row vector, we don't know if each element in that row vector
% corresponds to a different column or if the whole vector needs to be a
% property of each column.  Only the tech's MAT2CELL knows how to split the
% property appropriately amongst the columns.)
if ( (numel(f) == 1) && (numColumns(f) > 1) )
    out = getArbitraryLocalProp(cheb2quasi(f), prop);
    return
end

out = getArbitraryLocalPropQuasimatrix(f, prop);

% Try to convert from a cell array of cell arrays to a 2D cell array.
if ( doSimplify && canConvertCellOfCellsTo2DCellArray(out) )
    out = horzcat(out{:});

    % If our 2D cell array has only one cell (i.e., if f has only one column,
    % and that column has only fun), return the contents of that cell.
    if ( numel(out) == 1 )
        out = out{1};
    end
end

if ( f(1).isTransposed )
    out = out.';
end

end

function out = getSimpleNumericLocalPropArrayValued(f, prop)
%GETSIMPLENUMERICLOCALPROPARRAYVALUED   Get simple numeric FUN properties from
%   an array-valued CHEBFUN.
%
%   OUT = GETSIMPLENUMERICLOCALPROPARRAYVALUED(F, PROP) returns a cell array of
%   cell arrays OUT such that OUT{J}{K} is the value of the simple numeric
%   local property PROP (see GETSIMPLENUMERICLOCALPROP) for FUN J of column K
%   of the array-valued CHEBFUN F.

nCols = numColumns(f);
nFuns = numel(f.funs);

tmp = cell(nFuns, nCols);
for j = 1:nFuns
    funPropVal = get(f.funs{j}, prop);
    [funPropValRows, funPropValCols] = size(funPropVal);

    if ( funPropValCols == 1 )
        funPropVal = repmat(funPropVal, 1, nCols);
    elseif ( funPropValCols ~= nCols )
        error('CHEBFUN:get:notSimple', 'Property is not simple.');
    end

    tmp(j, :) = mat2cell(funPropVal, funPropValRows, ones(nCols, 1));
end

out = cell(1, nCols);
for k = 1:nCols
    out{k} = vertcat(tmp(:, k));
end

end

function out = getArbitraryLocalPropQuasimatrix(f, prop)
%GETARBITRARYLOCALPROPQUASIMATRIX   Get FUN properties from a quasimatrix.
%   OUT = GETARBITRARYLOCALPROPQUASIAMTRIX(F, PROP) returns a cell array of
%   cell arrays OUT such that OUT{J}{K} is the value of the local property PROP
%   for FUN J of column K of the quasimatrix F.

nCols = numColumns(f);
out = cell(1, nCols);
for k = 1:nCols
    nFuns = numel(f(k).funs);
    out{k} = cell(nFuns, 1);
    for j = 1:nFuns
        out{k}{j} = get(f(k).funs{j}, prop);
    end
end

end

function canConvert = canConvertCellOfCellsTo2DCellArray(out)
%CANCONVERTCELLOFCELLSTO2DCELLARRAY   Self-explanatory.
%   CANCONVERT = CANCONVERTCELLOFCELLTO2DCELLARRAY(OUT) function treturns TRUE
%   of the cell array of cell arrays OUT can be converted to a 2D cell array
%   (i.e., if OUT{J} has the same number of columns for each J) and FALSE
%   otherwise.

    canConvert = ~any(diff(cellfun(@(C) size(C, 1), out)));
end

function canConvert = canConvert2DCellArrayToNumeric(out)
%CANCONVERT2DCELLARRAYTONUMERIC   Self-explanatory.
%   CANCONVERT = CANCONVERT2DCELLARRAYTONUMERIC(OUT) function treturns TRUE of
%   the 2D cell array OUT can be converted to a numeric matrix.

[nFuns, nCols] = size(out);

% If each column has multiple funs, the property values for each fun must
% all have only one row.  If each column has only one fun, the property
% values must all have the same number of rows.
if ( nFuns > 1 )
    nRowsNeeded = 1;
else
    nRowsNeeded = size(out{1}, 1);
end

for ( j = 1:1:nFuns )
    for ( k = 1:1:nCols )
        if ( size(out{j, k}, 1) ~= nRowsNeeded )
            canConvert = false;
            return;
        end
    end
end

canConvert = true;

end
