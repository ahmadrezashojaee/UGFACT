function [model, states] = BioGeoModel_Parpool_Modified_Functionized(model, states, t, dt)
% BioGeoModel - Bio-Geochemical model - Considering no reaction when Sw = 0
%
% Syntax:
%       [model,states] = BioGeoModel(model,states,t)
%
% Description:
%       This function conduct a batch model in each grid after a timestep.
%       The batch model is developed in IPHREEQC.
%
% Inputs:
%       - model : The rock-fluid-grid model created in MRST as an input
%       - states: The matrix contains solution after each time step
%       - t     : Number of step
%       - dt    : length of time step
%
% Outputs:
%       - model : The updated model after batch calculations
%       - states: The updated solution matrix after batch calculations

state  = states{t,1};
G      = model.G;
rock   = model.rock;

% 1) Initialise biomass and geochem for this t
states = initBiomassAndGeochem(model, states, t);

% 2) Compute mineralogy moles in each cell
[states, mineralData] = computeMineralogy(model, states, t);

% 3) Compute phases, component masses, partial pressures etc.
phaseData = computePhaseMassesAndPressures(model, states, t, dt, mineralData);

% 4) Build PHREEQC strings
Strings = buildPhreeqcStrings(model, states, t, phaseData, mineralData);

% 5) Call PHREEQC in parallel and collect OUTPUT
OUTPUT  = runPhreeqcForAllCells(Strings, phaseData.Sw);

% 6) Postprocess OUTPUT and update state variables
[states, reactionRates] = updateStatesFromPhreeqc(model, states, t, OUTPUT, phaseData, mineralData);

% 7) Update MRST mass / Zi etc.
states = updateMassAndZi(model, states, t, reactionRates, phaseData);

end