classdef cheboppref < chebpref
% TODO: Add CHEBFUNPREF style documentation of the options available. To be written
% when we introduce nonlinear ODEs, since those involve a number of options.
    
% TODO: The relationship between CHEBOPPREF and CHEBFUNPREF needs some serious
% consideration.

    methods

        function outPref = cheboppref()           
            %outPref = outPref@chebfunpref;

            % Default new properties.
            outPref.prefList.domain = [-1 1];
            outPref.prefList.discretization = @colloc2;
            outPref.prefList.scale = NaN;
            outPref.prefList.dimensionValues = [32 64 128 256 512 724 1024 1448 2048];
            outPref.prefList.damped = 1;
            outPref.prefList.display = 'off';
            outPref.prefList.errTol = 1e-10;
            outPref.prefList.lambdaMin = 1e-6;
            outPref.prefList.maxIter = 25;
            outPref.prefList.plotting = 'off';
        end
        
    end

end
