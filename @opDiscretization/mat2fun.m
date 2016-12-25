function f = mat2fun(disc, values, cutoff)
%MAT2FUN    Convert a matrix of values to cell of CHEBFUNs.
%
%   Input data layout:
%
%     [   var1,piece1,col1       var1,piece1,col2      ...  ]
%     [      ...                     ...               
%     [   var1,pieceM,col1       var1,pieceM,col2      ...  ]
%     [      ...                     ...
%     [   varN,piece1,col1       varN,piece1,col2      ...  ]
%     [      ...                     ...  
%     [   varN,pieceM,col1       varN,pieceM,col2      ...  ]
%
%   The variables may be scalar valued. The input may be a cell, in which case it
%   will first be converted to the matrix form.
%
%   Output data layout is a cell vector. Each element is an array-valued
%   piecewise CHEBFUN, or a numeric row vector:
%
%     {   [ var1,col1       var1,col2      ...  ]   } 
%     {                                             }
%     {   [ var2,col1       var2,col2    ...    ]   }
%     {                                             }
%     {                   ...                       }
%     {                                             }
%     {   [ varN,col1       varN,col2      ...  ]   }

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Start by obtaining required information:
n = disc.dimension;                 % Discretization points of each function.
isFun = disc.source.isFunVariable;  % Which variables are functions?
numVar = length(isFun);             % Number of functions to return.

% If input is a cell, convert to matrix form:
if ( iscell(values) )
    values = cell2mat(values);
end

% Break each chebmatrix component into its own cell.:
componentLength = ones(1, numVar); % all scalars.
componentLength(isFun) = sum(n);   % total length of each function component.
values = mat2cell( values, componentLength, size(values,2) );

% Output is the same size cell, but with each function component entry converted
% from vector of values to chebfun:
f = cell(numVar, 1);
for j = 1:numVar
    if ( isFun(j) && nargin == 3 )
        f{j} = toFunctionOut(disc, values{j}, cutoff);
    elseif ( isFun(j) )
        f{j} = toFunctionOut(disc, values{j});
    else
        f{j} = values{j};
    end
end

end
