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

% Start by obtaining required information
isFun = disc.source.isFunVariable;  % Which variables are functions? 
numVar = length(isFun);             % Number of functions to return 

% Take discrete values and convert to CHEBFUN objects.
f = mat2fun(disc,values);

% Loop through the cell array C, and obtain the coefficients of all CHEBFUNS in
% C.
c = cell(numVar,1);
for i = 1:numVar
    if ( isFun(i) )
        c{i} = get(f{i}, 'coeffs').';
    else
        c{i} = f{i};
    end
end

end
