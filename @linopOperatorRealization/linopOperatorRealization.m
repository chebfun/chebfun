classdef (Abstract) linopOperatorRealization 

    methods (Abstract)
        
        % All operations are defined in terms of these three basics.
        C = plus(A,B)
        B = uminus(A)
        C = mtimes(A,B)

        % Given an (empty) instance of the realization object as the first argument,
        % the function should return the appropriate instantiation.
        D = diff(A,domain,m)      % differentiation
        C = cumsum(A,domain,m)    % indefinite integration
        I = eye(A,domain)       % identity
        Z = zeros(A,domain)     % zero
        F = diag(A,f)           % function multiplication
        % X = outer(A,f,g)
        
    end
          
end
  