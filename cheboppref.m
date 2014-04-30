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

       function display(pref)
       %DISPLAY   Display a CHEBFUNPREF object.
       %   DISPLAY(PREF) prints out a list of the preferences stored in the
       %   CHEBFUNPREF object PREF.

            % Compute the screen column in which pref values start.
            valueCol = 24; % length('    enableSingularityDetection:   ');

            % A subfunction to pad strings for formatting.
            function s = padString(s)
            %PADSTRING   Add whitespace to string for formatting.
                s = [s repmat(' ', 1, valueCol - length(s))];
            end

            % Print values of "known" preferences.
            prefList = pref.prefList;

            fprintf('cheboppref object with the following preferences:\n');
            fprintf([padString('    domain:') '[%g, %g]\n'], ...
                prefList.domain(1), prefList.domain(end));
            fprintf([padString('    discretization:') '%s\n'], ...
                func2str(prefList.discretization));
            fprintf([padString('    dimensionValues:') '%s\n'], ...
                num2str(prefList.dimensionValues));
            fprintf([padString('    damped:') '%d\n'], ...
                prefList.damped);
            fprintf([padString('    display:') '%s\n'], ...
                prefList.display);
            fprintf([padString('    errTol:') '%g\n'], ...
                prefList.errTol);
            fprintf([padString('    lambdaMin:') '%g\n'], ...
                prefList.lambdaMin);
            fprintf([padString('    maxIter:') '%d\n'], ...
                prefList.maxIter);
            fprintf([padString('    plotting:') '%s\n'], ...
                prefList.plotting);
        end

    end

end
