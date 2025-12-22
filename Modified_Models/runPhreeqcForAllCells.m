function OUTPUT = runPhreeqcForAllCells(Strings, activeCells)
%RUNPHREEQCFORALLCELLS
% Run PHREEQC for active cells only.
% Prints cell numbers that did not converge.

    nc     = numel(Strings);
    OUTPUT = cell(nc, 1);

    % Database
    dbPath = 'database\PHREEQC_Modified.dat';

    % Logical flag: true if this cell failed
    errorFlag = false(nc, 1);

    for i = 1:nc
        OUTPUT{i} = [];  % default

        if activeCells(i) && ~isempty(Strings{i})
            hasError = false;
            try
                iph = actxserver('IPhreeqcCOM.Object');
                iph.LoadDatabase(dbPath);

                rc = iph.RunString(Strings{i});   % 0 = success

                if rc ~= 0
                    % PHREEQC reported an error
                    hasError = true;
                else
                    out = iph.GetSelectedOutputArray;
                    if isempty(out)
                        % Empty output → treat as error
                        hasError = true;
                    else
                        OUTPUT{i} = out;
                    end
                end
            catch
                % MATLAB / COM error
                hasError = true;
            end

            errorFlag(i) = hasError;
        else
            errorFlag(i) = false;
        end
    end

    % ------------ Print only failing cell numbers -------------
    failedIdx = find(errorFlag);

    if ~isempty(failedIdx)
        fprintf('PHREEQC failed in %d cells:\n', numel(failedIdx));
        fprintf('Cells: %s\n', num2str(failedIdx(:).'));
    else
        fprintf('PHREEQC: all active cells converged.\n');
    end
end

