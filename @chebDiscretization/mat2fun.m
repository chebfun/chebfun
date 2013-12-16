function f = mat2fun(disc,values)
% Convert a matrix of values to cell of chebfuns.
%
% Input data layout:
%
%     [   var1,piece1,col1       var1,piece1,col2      ...  ]
%     [      ...                     ...               
%     [   var1,pieceM,col1       var1,pieceM,col2      ...  ]
%     [      ...                     ...
%     [   varN,piece1,col1       varN,piece1,col2      ...  ]
%     [      ...                     ...  
%     [   varN,pieceM,col1       varN,pieceM,col2      ...  ]
%
% The variables may be scalar valued.
%
% Output data layout is a cell vector. Each element is an array-valued piecewise
% chebfun, or a numeric row vector:
%
%     {   [ var1,col1       var1,col2      ...  ]   } 
%     {                                             }
%     {   [ var2,col1       var2,col2    ...    ]   }
%     {                                             }
%     {                   ...                       }
%     {                                             }
%     {   [ varN,col1       varN,col2      ...  ]   }
%

n = disc.dimension;
isFun = disc.source.isFunVariable; 
numvar = length(isFun);

m = ones(1,numvar);
m(isFun) = sum(n);
values = mat2cell( values, m, size(values,2) );

f = cell(numvar,1);
for j = 1:numvar
    if ( isFun(j) )
        f{j} = toFunction(disc,values{j});
    else
        f{j} = values{j};
    end
end

end
