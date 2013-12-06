function c = mat2poly(disc,values)
% Convert a matrix of values to cell of Chebyshev polynomial coefficients.
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
%     {  {  [ var1,piece1,col1   var1,piece1,col2  ...  ]  }  } 
%     {  {                     ...                         }  } 
%     {  {  [ var1,pieceM,col1   var1,pieceM,col2  ...  ]  }  } 
%     {                      ...                              }
%     {  {  [ varN,piece1,col1   varN,piece1,col2  ...  ]  }  } 
%     {  {                     ...                         }  } 
%     {  {  [ varN,pieceM,col1   varN,pieceM,col2  ...  ]  }  } 
%
% That is, out{i}{j} is an array of polynomial coefficients for the ith variable
% on the jth interval, with rows giving decreasing degree. If the ith variable
% is a scalar, then out{i} is just a row vector of values. 

isFun = disc.source.isFunVariable; 
numvar = length(isFun);

f = mat2fun(disc,values);
c = cell(numvar,1);
for i = 1:numvar
    if ( isFun(i) )
        for j = 1:disc.numIntervals
            c{i}{j} = chebpoly(f{i},j).';
        end
    else
        c{i} = f{i};
    end
end

end
