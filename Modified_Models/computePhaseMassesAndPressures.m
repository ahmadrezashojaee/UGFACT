function phaseData = computePhaseMassesAndPressures(model, states, t, dt, mineralData) %#ok<INUSD>
%COMPUTEPHASEMASSESANDPRESSURES
% Compute water/gas mixture properties, component masses and moles,
% partial pressures, and cleaned aqueous concentrations for step t.
%
% Inputs:
%   model       - MRST model
%   states      - state cell array
%   t           - current time step index
%   dt          - time step length (seconds)
%   mineralData - kept for interface consistency (not used here)
%
% Output:
%   phaseData   - struct with fields used later by PHREEQC string builder
%                 and post-processing (masses, moles, Sw/Sg, pH, K, Na, ...)

state = states{t,1};
G     = model.G;
rock  = model.rock;
nc    = G.cells.num;

% ---- Mixture densities ----
Water_Density = state.PVTProps.Density{1,1};  % [kg/m^3]
Gas_Density   = state.PVTProps.Density{1,2};  % [kg/m^3]

% Overall compositions
X = state.x;           % aqueous overall mole fractions
Y = state.y;           % gas overall mole fractions
Z = state.components;  % overall cell composition z_i

% Molecular weights [kg/mol]
MW_H2O = 0.018015268;
MW_H2  = 0.00201588;
MW_CO2 = 0.0440098;
MW_CH4 = 0.0160428;
MW_H2S = 0.03408088;
MW_N2  = 0.02801348;

% Mixture molecular weights
MW_Water = X(:,1).*MW_H2O + X(:,2).*MW_H2  + X(:,3).*MW_CO2 + ...
    X(:,4).*MW_CH4 + X(:,5).*MW_H2S + X(:,6).*MW_N2;

MW_Gas   = Y(:,1).*MW_H2O + Y(:,2).*MW_H2  + Y(:,3).*MW_CO2 + ...
    Y(:,4).*MW_CH4 + Y(:,5).*MW_H2S + Y(:,6).*MW_N2;

rho_molar_Water = Water_Density ./ MW_Water;  % [mol/m^3]
rho_molar_Gas   = Gas_Density   ./ MW_Gas;    % [mol/m^3]

% ---- Saturations ----
Sw = state.s(:,1);
Sg = state.s(:,2);

% ---- Water mass per cell (kg) ----
mass_H2O_Profile = Sw .* G.cells.volumes .* rock.poro .* ...
    rho_molar_Water .* X(:,1) .* MW_H2O;   % [kg]

% ---- Component masses from MRST (already in kg) ----
mass_H2_Profile  = state.FlowProps.ComponentTotalMass{2,1};
mass_CO2_Profile = state.FlowProps.ComponentTotalMass{3,1};
mass_CH4_Profile = state.FlowProps.ComponentTotalMass{4,1};
mass_H2S_Profile = state.FlowProps.ComponentTotalMass{5,1};
mass_N2_Profile  = state.FlowProps.ComponentTotalMass{6,1};

% Component overall fractions (for tiny check)
H2component_Profile  = Z(:,2);
CO2component_Profile = Z(:,3);
CH4component_Profile = Z(:,4);
H2Scomponent_Profile = Z(:,5);
N2component_Profile  = Z(:,6);

% ---- Aqueous chemistry and tolerances ----
pH    = state.Solution.pH;
Units = state.Solution.Unit;

tol = 1e-6;

K       = state.Solution.K;
Na      = state.Solution.Na;
Mg      = state.Solution.Mg;
Ca      = state.Solution.Ca;
Cl      = state.Solution.Cl;
C4      = state.Solution.C4;
S6      = state.Solution.S6;
S2      = state.Solution.S2;
Fe3     = state.Solution.Fe3;
Fe2     = state.Solution.Fe2;
Si      = state.Solution.Si;
Acetate = state.Solution.Acetate;

% Zero out tiny negative / tiny values
K(K < tol)         = 0;
Na(Na < tol)       = 0;
Mg(Mg < tol)       = 0;
Ca(Ca < tol)       = 0;
Cl(Cl < tol)       = 0;
C4(C4 < tol)       = 0;
S6(S6 < tol)       = 0;
S2(S2 < tol)       = 0;
Fe3(Fe3 < tol)     = 0;
Fe2(Fe2 < tol)     = 0;
Si(Si < tol)       = 0;
Acetate(Acetate < tol) = 0;

% Biomass (carried from initBiomassAndGeochem)
SRB_Biomass = state.Solution.SRB_Biomass;
MET_Biomass = state.Solution.MET_Biomass;
ACE_Biomass = state.Solution.ACE_Biomass;
FRB_Biomass = state.Solution.FRB_Biomass;

% Gas constant in atm·L/(mol·K)
R = 0.082057366080960;

% ---- Clean tiny gas fractions (your original logic) ----
idxTinyH2   = H2component_Profile  <= 1e-8;
idxTinyCO2  = CO2component_Profile <= 1e-8;
idxTinyCH4  = CH4component_Profile <= 1e-8;
idxTinyH2S  = H2Scomponent_Profile <= 1e-8;
idxTinyN2   = N2component_Profile  <= 1e-8;

mass_H2_Profile(idxTinyH2)   = 0;
mass_CO2_Profile(idxTinyCO2) = 0;
mass_CH4_Profile(idxTinyCH4) = 0;
mass_H2S_Profile(idxTinyH2S) = 0;
mass_N2_Profile(idxTinyN2)   = 0;

% This index is used later when you decide to keep old pressure
idxKeepOldPressure = idxTinyN2;

% ---- Moles of water and gas components ----
mol_H2O_Profile = mass_H2O_Profile ./ MW_H2O;
mol_H2_Profile  = mass_H2_Profile  ./ MW_H2;
mol_CO2_Profile = mass_CO2_Profile ./ MW_CO2;
mol_CH4_Profile = mass_CH4_Profile ./ MW_CH4;
mol_H2S_Profile = mass_H2S_Profile ./ MW_H2S;
mol_N2_Profile  = mass_N2_Profile  ./ MW_N2;

% ---- Gas volume (L) ----
Vg = G.cells.volumes .* Sg .* 1000 .* rock.poro;  % [L]
idxZeroVg = (Vg == 0);
Vg(idxZeroVg) = 1e-6;    % small epsilon

% ---- Partial pressures (atm) ----
T = state.T;   % [K]
pH2  = state.Z_V .* mol_H2_Profile  .* R .* T ./ Vg;
pCO2 = state.Z_V .* mol_CO2_Profile .* R .* T ./ Vg;
pCH4 = state.Z_V .* mol_CH4_Profile .* R .* T ./ Vg;
pH2S = state.Z_V .* mol_H2S_Profile .* R .* T ./ Vg;
pN2  = state.Z_V .* mol_N2_Profile  .* R .* T ./ Vg;

% ---- Time and PHREEQC integration parameters ----
tot_time = dt;     % [s]
if isfield(model,'BioGeoSteps') && model.BioGeoSteps>=1
    steps = model.BioGeoSteps;
else
    steps    = 3;      % number of internal PHREEQC time steps
end

% ---- Pressure for PHREEQC (atm) and temperature (K) ----
pressure_atm = state.pressure ./ 101325;  % Pa → atm
temp         = T;

% ---- Active mask for PHREEQC based on H2 and saturation ----
zH2    = Z(:,2);      % overall H2 mole fraction
tolH2  = 1e-3;

if (t == 1 && model.state0.initGeoChem)
    % At the first time step, run PHREEQC for all cells
    activeCells = true(model.G.cells.num, 1);
else
    % Normal rule: must have both water AND hydrogen
    hasWater = Sw >= 1e-3;
    hasH2    = zH2 > tolH2;

    activeCells = hasWater & hasH2;
end


% ---- Pack everything into a struct ----
phaseData = struct( ...
    'Water_Density',      Water_Density, ...
    'Gas_Density',        Gas_Density, ...
    'X',                  X, ...
    'Y',                  Y, ...
    'Z',                  Z, ...
    'MW_H2O',             MW_H2O, ...
    'MW_H2',              MW_H2, ...
    'MW_CO2',             MW_CO2, ...
    'MW_CH4',             MW_CH4, ...
    'MW_H2S',             MW_H2S, ...
    'MW_N2',              MW_N2, ...
    'MW_Water',           MW_Water, ...
    'MW_Gas',             MW_Gas, ...
    'rho_molar_Water',    rho_molar_Water, ...
    'rho_molar_Gas',      rho_molar_Gas, ...
    'mass_H2O_Profile',   mass_H2O_Profile, ...
    'mass_H2_Profile',    mass_H2_Profile, ...
    'mass_CO2_Profile',   mass_CO2_Profile, ...
    'mass_CH4_Profile',   mass_CH4_Profile, ...
    'mass_H2S_Profile',   mass_H2S_Profile, ...
    'mass_N2_Profile',    mass_N2_Profile, ...
    'mol_H2O_Profile',    mol_H2O_Profile, ...
    'mol_H2_Profile',     mol_H2_Profile, ...
    'mol_CO2_Profile',    mol_CO2_Profile, ...
    'mol_CH4_Profile',    mol_CH4_Profile, ...
    'mol_H2S_Profile',    mol_H2S_Profile, ...
    'mol_N2_Profile',     mol_N2_Profile, ...
    'Sw',                 Sw, ...
    'Sg',                 Sg, ...
    'pH',                 pH, ...
    'Units',              Units, ...
    'K',                  K, ...
    'Na',                 Na, ...
    'Mg',                 Mg, ...
    'Ca',                 Ca, ...
    'Cl',                 Cl, ...
    'C4',                 C4, ...
    'S6',                 S6, ...
    'S2',                 S2, ...
    'Fe3',                Fe3, ...
    'Fe2',                Fe2, ...
    'Si',                 Si, ...
    'Acetate',            Acetate, ...
    'SRB_Biomass',        SRB_Biomass, ...
    'MET_Biomass',        MET_Biomass, ...
    'ACE_Biomass',        ACE_Biomass, ...
    'FRB_Biomass',        FRB_Biomass, ...
    'R',                  R, ...
    'Vg',                 Vg, ...
    'pH2',                pH2, ...
    'pCO2',               pCO2, ...
    'pCH4',               pCH4, ...
    'pH2S',               pH2S, ...
    'pN2',                pN2, ...
    'tot_time',           tot_time, ...
    'steps',              steps, ...
    'pressure_atm',       pressure_atm, ...
    'temp',               temp, ...
    'idxKeepOldPressure', idxKeepOldPressure, ...
    'zH2',                zH2, ...
    'tolH2',              tolH2, ...
    'activeCells',        activeCells, ...
    'nc',                 nc ...
    );
end
