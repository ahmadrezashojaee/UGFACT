function states = updateMassAndZi_Parpool(model, states, t, reactionRates, phaseData) %#ok<INUSD>
%UPDATEMASSANDZI
% Final post-processing step after PHREEQC:
% assumes FlowProps.ComponentTotalMass has already been updated
% (for active cells) and simply recomputes Zi etc. via mass2Zi.
%
% Inputs:
%   model        - MRST model
%   states       - state cell array
%   t            - current time step index
%   reactionRates- struct from updateStatesFromPhreeqc (kept for interface)
%   phaseData    - struct from computePhaseMassesAndPressures (not used here)
%
% Output:
%   states       - state cell array with updated Zi in states{t,1}

    states = mass2Zi_Parpool(states, model, t);

end
