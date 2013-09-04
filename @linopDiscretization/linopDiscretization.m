classdef (Abstract) linopDiscretization < linopOperatorRealization & linopFunctionalRealization
  
    properties (Abstract)
        size
    end
    
    methods (Abstract,Static)
        B = resize(A,m)
        isDone = convergeTest(v)
    end
    
end
        