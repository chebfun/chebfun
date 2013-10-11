classdef (Abstract) functionalBlockRealization 

    methods (Abstract)

        % All operations are defined in terms of these three basics.
        C = plus(A,B)
        B = uminus(A)
        C = mtimes(A,B)

        % Given an (empty) instance of the realization object as the first argument,
        % the function should return the appropriate instantiation.
        Z = zero(A,domain)
        S = sum(A,domain)         % definite integration     
        E = feval(A,domain,loc)   % point evaluation
        F = inner(A,f)            % function inner product                
        
    end
         
end
  