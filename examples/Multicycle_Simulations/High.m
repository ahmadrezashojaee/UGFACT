
%% Set up problem
% Define grid
clear;clc;close all
mrstModule add UGFACT2 ad-core ad-props mrst-gui ad-blackoil deckformat
%% Gridding
fn = 'grid.GRDECL';
G = readGRDECL(fn);
G = processGRDECL(G);
G = computeGeometry(G);
NX = G.cartDims(1);
NY = G.cartDims(2);
NZ = G.cartDims(3);
Lx = max(G.nodes.coords(:,1)) - min(G.nodes.coords(:,1));
Ly = max(G.nodes.coords(:,2)) - min(G.nodes.coords(:,2));
Lz = max(G.nodes.coords(:,3)) - min(G.nodes.coords(:,3));
%% Rock model
%% Top Section
% Correlation lengths (m)
corrX1=300; corrY1=corrX1; corrZ1=1;

% Porosity stats
poro_median1=0.10; poro_std1=0.01; poro_min1=0.05;

% Permeability stats (mD)
perm_median1=65; perm_std1=10;

rho1=0.90;    % poro-perm correlation strength
seed1=14000;  % reproducibility

[poro1, perm1] = generateReservoirPropertiesFFTMA(NX,NY,NZ/2,Lx,Ly,Lz/2, ...
    corrX1,corrY1,corrZ1, ...
    poro_median1,poro_std1,poro_min1, ...
    perm_median1,perm_std1,rho1,seed1);
%% Bottom Section
% Correlation lengths (m)
corrX2=300; corrY2=corrX2; corrZ2=3;

% Porosity stats
poro_median2=0.14; poro_std2=0.02; poro_min2=0.08;

% Permeability stats (mD)
perm_median2=85; perm_std2=10;

rho2=0.95;    % poro-perm correlation strength
seed2=14000;  % reproducibility
[poro2, perm2] = generateReservoirPropertiesFFTMA(NX,NY,NZ/2,Lx,Ly,Lz/2, ...
    corrX2,corrY2,corrZ2, ...
    poro_median2,poro_std2,poro_min2, ...
    perm_median2,perm_std2,rho2,seed2);
%% Merge them
poro = [poro1;poro2];
perm = [perm1;perm2];
% Cross-plot
figure;
scatter(poro,perm,10,'filled');
xlabel('Porosity [-]'); ylabel('Permeability [mD]');
title('Porosity–Permeability Cross-plot');
grid on;


% Plot permeability cross-section
xc = reshape(G.cells.centroids(:,1), G.cartDims);   % X of each cell
zc = reshape(G.cells.centroids(:,3), G.cartDims);   % Z of each cell

X = squeeze(xc(:,1,:))';    % [NZ, NX]
Z = squeeze(zc(:,1,:))';    % [NZ, NX]
perm2D = reshape(perm,[NX,NZ])';
poro2D = reshape(poro,[NX,NZ])';
figure;
% --- Permeability ---
subplot(1,2,2)
surf(X, Z, perm2D, 'EdgeColor', 'none');
view(2); shading interp; colormap(turbo); colorbar;
set(gca, 'YDir', 'reverse');
xlabel('X (m)'); ylabel('Z (m)');
title('Permeability [mD]');
axis tight                % <--- remove axis equal
% --- Porosity ---
subplot(1,2,1)
surf(X, Z, poro2D, 'EdgeColor', 'none');
view(2); shading interp; colormap(turbo); colorbar;
set(gca, 'YDir', 'reverse');
xlabel('X (m)'); ylabel('Z (m)');
title('Porosity [fraction]');
axis tight                % <--- remove axis equal

% initial_k = 1000*milli*darcy;
% initial_phi = 0.2;
pv = sum(G.cells.volumes.*poro);
rock = makeRock(G, [perm perm 0.1*perm]*milli*darcy, poro); % Making a rock model based on Grid, perm and poro
%% Fluid model
f = initSimpleADIFluid('phases', 'og', 'blackoil', false);
[WaterTable, GasTable] = coreyBrooksTables(0.24, 0.05, 1, 1, 3, 2, 10);
f.krO = @(sw)interpTable(WaterTable(:,1),WaterTable(:,2),sw);
f.krG = @(sg)interpTable(GasTable(:,1),GasTable(:,2),sg);
mixture = TableCompositionalMixture({'Water','Hydrogen','CarbonDioxide','Methane','HydrogenSulfide','Nitrogen'},...
    {'H2O','H2','CO2','CH4','H2S','N2'});

bic = zeros(6,6);
%bicCO2H2O = -0.31092*(1+0.15587*0.001^0.7505) + 0.23580*(1+0.17837*0.001^0.979) - 21.2566*exp(-6.7222*1.2-0.001);
bic(1,3) = -0.04;
bic(3,1) = bic(1,3);
bic(1,4) = -0.12;
bic(4,1) = bic(1,4);
bic(1,2) = -0.58;
bic(2,1) = bic(1,2);
mixture = mixture.setBinaryInteraction(bic);
arg = {G, rock, f, ... % Standard arguments
    mixture,... % Compositional mixture
    'water', false, 'oil', true, 'gas', true,... % Water-Gas system
    'liquidPhase', 'O', 'vaporPhase', 'G'}; % Water=liquid, gas=vapor
% Construct models for both formulations. Same input arguments

%% Defining the model based on the Grid, Rock, fluid, and mixture.
model = GenericOverallCompositionModel(arg{:}); % Overall mole fractions model

%model = GenericNaturalVariablesModel(arg{:}); % Natural variables
model.EOSModel.PropertyModel.volumeShift = [0, 0, 0, 0, 0, 0]; % Volume shift
gravity reset on
%% Initial state
p = 150*barsa; T = 273.15 + 92; s = []; z = [0.7, 0, 0, 0.3, 0, 0]; % p, T, s, z

state0   = initCompositionalState(model, p, T, s, z); % Initialize state
%% Input data for biochemical and geochemical reactions
model.BioGeo = true;
% model.parpool = false;
% model.parpoolCores  = 5;
% model.BioGeoSteps  = 5;
model.Solution.pH    = repmat(6.9,G.cells.num,1);
model.Solution.Unit  = 'mol/kgw';
%Ion concentration in the given unit
model.Solution.Ca    = repmat(2.368e-02,G.cells.num,1);
model.Solution.Cl    = repmat(1.133e-01,G.cells.num,1);
model.Solution.Na    = repmat(1.039e-01,G.cells.num,1);
model.Solution.K     = repmat(4.998e-03,G.cells.num,1); %#ok<RPMT0>
model.Solution.S6    = repmat(2.013e-03,G.cells.num,1);
model.Solution.C4    = repmat(3.338e-03,G.cells.num,1);
model.Solution.Mg    = repmat(1.826e-02,G.cells.num,1);
model.Solution.Fe3   = repmat(0,G.cells.num,1); %#ok<REPMAT>
model.Solution.Fe2   = repmat(1.694e-13,G.cells.num,1); %#ok<REPMAT>
model.Solution.Acetate   = repmat(0,G.cells.num,1); %#ok<REPMAT>
model.Solution.S2    = repmat(1.201e-13,G.cells.num,1); %#ok<REPMAT>;
model.Solution.Si    = repmat(9.723e-05,G.cells.num,1); %#ok<REPMAT>;

model.Mineralogy.Calcite      = repmat(0,G.cells.num,1);  %#ok<REPMAT> %Percent
model.Mineralogy.Dolomite     = repmat(2.42,G.cells.num,1);  %#ok<REPMAT> % %Percent
model.Mineralogy.Dolomite(251:500) = 0.1; %Top Section of reservoir has Dolomite
model.Mineralogy.Quartz       = repmat(97.53,G.cells.num,1); %#ok<REPMAT> % %Percent
model.Mineralogy.Quartz(251:500) = 99.85;
model.Mineralogy.Anhydrite    = repmat(0.05,G.cells.num,1); %#ok<REPMAT> %Percent
model.Mineralogy.Goethite     = repmat(0,G.cells.num,1); %#ok<REPMAT> %Percent
model.Mineralogy.Brucite      = repmat(0,G.cells.num,1); %#ok<REPMAT> %Percent
model.Mineralogy.Portlandite  = repmat(0,G.cells.num,1); %#ok<REPMAT> %Percent
model.Mineralogy.Pyrite       = repmat(0,G.cells.num,1); %#ok<REPMAT> %Percent
model.Mineralogy.Gypsum       = repmat(0,G.cells.num,1); %#ok<REPMAT> %Percent

%% Kinetic Data
%Methanogenesis
model.Kinetic.mu_MET   = 4.1;%1.109;  0.3 4.1      %per day
model.Kinetic.b_MET    = 0.01 * model.Kinetic.mu_MET; %Decay Coefficent
model.Kinetic.Y_MET    = 0.04;       %Growth Yield
model.Kinetic.K_Dmet   = 9e-6;      %Electron Donor Half Saturation Constant
model.Kinetic.K_Amet   = 230e-6;     %Electron Acceptor Half Saturation Constant
model.Kinetic.MW_MET   = 24.6;         %Microbe Molecular Weight g/mol
model.Kinetic.m_MET    = 1e-14;        %Biomass Mass gram
model.Kinetic.N0_MET   = 1e6*1000;           %number of biomass per 1 kg of water
model.Kinetic.M0_MET   = model.Kinetic.N0_MET * model.Kinetic.m_MET * 1/model.Kinetic.MW_MET; %Initial mole of Biomass
model.Kinetic.Nmax_MET = 1e10*1000;       %Maximum biomass concnetration
model.Kinetic.Mmax_MET = model.Kinetic.Nmax_MET * model.Kinetic.m_MET * 1/model.Kinetic.MW_MET; %Maximum biomass mole per 1 kg of water
model.Kinetic.Mmin_MET = model.Kinetic.N0_MET * model.Kinetic.m_MET / model.Kinetic.MW_MET; %Minimum biomass mole per 1 kg of water
%Sulfate reduction
model.Kinetic.mu_SRB   = 5.5;%5.5;%1.048; 0.2       %per day
model.Kinetic.b_SRB    = 0.01 * model.Kinetic.mu_SRB; %Decay Coefficent
model.Kinetic.Y_SRB    = 0.11;       %Growth Yield
model.Kinetic.K_Dsrb   = 2.9e-6;       %Electron Donor Half Saturation Constant
model.Kinetic.K_Asrb   = 2751.5e-6;    %Electron Acceptor Half Saturation Constant
model.Kinetic.MW_SRB   = 24.6;         %Microbe Molecular Weight g/mol
model.Kinetic.m_SRB    = 1e-14;        %Biomass Mass gram
model.Kinetic.N0_SRB   = 1e6*1000;           %numnber of biomass per 1 kg of water
model.Kinetic.M0_SRB   = model.Kinetic.N0_SRB * model.Kinetic.m_SRB * 1/model.Kinetic.MW_SRB; %Initial mole of Biomass per 1kg of water
model.Kinetic.Nmax_SRB = 1e10*1000;       %Maximum biomass concnetration
model.Kinetic.Mmax_SRB = model.Kinetic.Nmax_SRB * model.Kinetic.m_SRB * 1/model.Kinetic.MW_SRB; %Maximum biomass mole per 1 kg of water
model.Kinetic.Mmin_SRB = model.Kinetic.N0_SRB * model.Kinetic.m_SRB / model.Kinetic.MW_SRB; %Minimum biomass mole per 1 kg of water
%Acetogenesis
model.Kinetic.mu_ACE   = 0;%1.9;%0.872; 0.4       %per day
model.Kinetic.b_ACE    = 0.01 * model.Kinetic.mu_ACE; %Decay Coefficent
model.Kinetic.Y_ACE    = 0.07;       %Growth Yield
model.Kinetic.K_Dace   = 2.5e-6;       %Electron Donor Half Saturation Constant
model.Kinetic.K_Aace   = 115.5e-6;     %Electron Acceptor Half Saturation Constant
model.Kinetic.MW_ACE   = 24.6;         %Microbe Molecular Weight g/mol
model.Kinetic.m_ACE    = 1e-14;        %Biomass Mass gram
model.Kinetic.N0_ACE   = 1e6*1000;           %numnber of biomass per 1 kg of water
model.Kinetic.M0_ACE   = model.Kinetic.N0_ACE * model.Kinetic.m_ACE * 1/model.Kinetic.MW_ACE; %Initial mole of Biomass
model.Kinetic.Nmax_ACE = 1e10*1000;       %Maximum biomass concnetration
model.Kinetic.Mmax_ACE = model.Kinetic.Nmax_ACE * model.Kinetic.m_ACE * 1/model.Kinetic.MW_ACE; %Maximum biomass mole per 1 kg of water
model.Kinetic.Mmin_ACE = model.Kinetic.N0_ACE * model.Kinetic.m_ACE / model.Kinetic.MW_ACE; %Minimum biomass mole per 1 kg of water
%Iron reduction
model.Kinetic.mu_FRB   = 0;%1.5;        %per day
model.Kinetic.b_FRB    = 0.01 * model.Kinetic.mu_FRB; %Decay Coefficent
model.Kinetic.Y_FRB    = 0.14*4;       %Growth Yield
model.Kinetic.K_Dfrb   = 1e-6;       %Electron Donor Half Saturation Constant
model.Kinetic.K_Afrb   = 1e-12;     %#ok<NASGU> %Electron Acceptor Half Saturation Constant
model.Kinetic.MW_FRB   = 24.6;         %Microbe Molecular Weight g/mol
model.Kinetic.m_FRB    = 1e-14;        %Biomass Mass gram
model.Kinetic.N0_FRB   = 1e6*1000;           %numnber of biomass per 1 kg of water
model.Kinetic.M0_FRB   = model.Kinetic.N0_FRB * model.Kinetic.m_FRB * 1/model.Kinetic.MW_FRB; %Initial mole of Biomass
model.Kinetic.Nmax_FRB = 1e10*1000;       %Maximum biomass concnetration
model.Kinetic.Mmax_FRB = model.Kinetic.Nmax_FRB * model.Kinetic.m_FRB * 1/model.Kinetic.MW_FRB; %Maximum biomass mole per 1 kg of water
model.Kinetic.Mmin_FRB = model.Kinetic.N0_FRB * model.Kinetic.m_FRB / model.Kinetic.MW_FRB; %Minimum biomass mole per 1 kg of water
model.state0.initGeoChem = true;
if (model.BioGeo && model.state0.initGeoChem)
    model = GeochemistryInitializer(model,state0);
end
% stand alone flash - this section is not necessary. Here is to find Z_V
eos = EquationOfStateModel([], mixture, 'Peng-Robinson');
[L, x, y, Z_L, Z_V, rho_L, rho_G] = standaloneFlash(p, T, [0,1,0,0,0,0], eos);

%% Storage Scenario


Bg = 101325/298.15*Z_V*T/p;
rate = 0.002*pv*meter^3/day/Bg; % Surface Rate
% Parameters
% Well controls
k_range = 1:3;
j = 1; i = 25;   % first layer
cells = sub2ind([NX NY NZ], i*ones(size(k_range)), j*ones(size(k_range)), k_range);
% Build schedule and wells
injType='grat'; injVal=rate;
prodType='grat'; prodVal=-rate;

dt = 5*day;   % user-chosen timestep length

[schedule,W,stepDates] = makeAnnualStorageSchedule525(G, rock, cells, 2030, 2060, ...
                                       injType, injVal, prodType, prodVal, dt);

s  = EOSSeparator('pressure', 1*atm, 'T', 298.15); % Set conditions for surface
sg = SeparatorGroup(s);                         % Group = single separator
sg.mode = 'moles';                              % Use mole mode
model.FacilityModel.SeparatorGroup = sg;    % Connect to reservoir model
%% Running the simulation with a visualization tool
tic
[wellSol, states, report] = simulateScheduleAD_Modified(state0, model, schedule);
toc
