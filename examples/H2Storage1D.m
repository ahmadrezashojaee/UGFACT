%% Set up problem
% Define grid
clear;clc;close all
mrstModule add UGFACT ad-core ad-props mrst-gui ad-blackoil

%% Gridding
[NX,NY,NZ]=deal(50,1,1); % Number of cells in x,y,z direction
[dx,dy,dz]=deal(50,1,1); % Reservoir dimenstions in x,y,z direction 
G=computeGeometry(cartGrid([NX,NY,NZ],[dx,dy,dz])); % Making a cartesian grid
%% Rock model
initial_k = 100*milli*darcy;
initial_phi = 0.2;
rock = makeRock(G, initial_k, initial_phi); % Making a rock model based on Grid, perm and poro
%plotCellData(G,log10(convertTo(rock.perm(:,1),milli*darcy)))
%view(3)
pv = sum(G.cells.volumes.*rock.poro); % Calculating the pore volume of each cell
%% Fluid model
f = initSimpleADIFluid('phases', 'og', 'blackoil', false);
% SW_table = [-0.1 0:0.05:1 2];
% pc_table = [500 500 300 200 150 10 5 3.5 3 2 1.9 1.8 1.7 1.6 1.5 1.4 1.3 1.2 0.2 -0.1 -0.2 -0.3 -0.3];
% f.pcGW = @(sw)interpTable(SW_table,pc_table,sw); 
Sw_table = [0 0.16,0.20,0.24,0.28,0.32,0.36,0.4,0.44,0.48,0.52,0.56,0.6,0.64,0.68,0.72,0.76,0.8,0.84,0.88,0.92,0.96,0.999];
krW_table = [0 0 0.002,0.01,0.02,0.033,0.049,0.066,0.09,0.119,0.15,0.186,0.227,0.277,0.33,0.39,0.462,0.54,0.62,0.71,0.8,0.9,1];
f.krO = @(sw)interpTable(Sw_table,krW_table,sw);
Sg_table =[0.05,0.08,0.12,0.16,0.2,0.24,0.28,0.32,0.36,0.4,0.44,0.48,0.52,0.56,0.6,0.64,0.68,0.72,0.76,0.8,0.84 1];
krG_table =[0,0.013,0.026,0.04,0.058,0.078,0.1,0.126,0.156,0.187,0.222,0.26,0.3,0.348,0.4,0.45,0.505,0.562,0.62,0.68,0.74 1];
f.krG = @(sg)interpTable(Sg_table,krG_table,sg);
mixture = TableCompositionalMixture({'Water','Hydrogen','CarbonDioxide','Methane','HydrogenSulfide','Nitrogen'},...
                                    {'H2O','H2','CO2','CH4','H2S','N2'});

bic = zeros(6,6);
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
p = 150*barsa; T = 273.15 + 60; s = []; z = [0.9, 0, 0, 0.1, 0, 0]; % p, T, s, z

state0   = initCompositionalState(model, p, T, s, z); % Initialize state
%% Input data for biochemical and geochemical reactions
model.BioGeo = true;
model.parpool = false;   % Use parallel calculations
model.parpoolCores  = 5; % Number of cores
model.BioGeoSteps  = 5;  % Number of sub steps for bio-geochemical simulations 
model.Solution.pH    = repmat(6.24,G.cells.num,1);
model.Solution.Unit  = 'mol/kgw'; %ppm, mmol/kgw, 
%Ion concentration in the given unit
model.Solution.Ca    = repmat(2.857e-01,G.cells.num,1);
model.Solution.Cl    = repmat(3.655,G.cells.num,1);
model.Solution.Na    = repmat(2.865,G.cells.num,1);
model.Solution.K     = repmat(0,G.cells.num,1); %#ok<RPMT0>
model.Solution.S6    = repmat(4.664e-3,G.cells.num,1);
model.Solution.C4    = repmat(1.370e-03,G.cells.num,1); 
model.Solution.Mg    = repmat(1.144e-01,G.cells.num,1);
model.Solution.Fe3   = repmat(0,G.cells.num,1); %#ok<REPMAT>
model.Solution.Fe2   = repmat(0,G.cells.num,1); %#ok<REPMAT>
model.Solution.Acetate   = repmat(0,G.cells.num,1); %#ok<REPMAT>
model.Solution.S2    = repmat(0,G.cells.num,1); %#ok<REPMAT>;
model.Solution.Si    = repmat(9.723e-05,G.cells.num,1); %#ok<REPMAT>;

model.Mineralogy.Calcite      = repmat(0,G.cells.num,1); %#ok<REPMAT> % Weight Percentage
model.Mineralogy.Dolomite     = repmat(2,G.cells.num,1); %#ok<REPMAT> % Weight Percentage
model.Mineralogy.Quartz       = repmat(98,G.cells.num,1); %#ok<REPMAT> % Weight Percentage
model.Mineralogy.Anhydrite    = repmat(0,G.cells.num,1); %#ok<REPMAT> % Weight Percentage
model.Mineralogy.Goethite     = repmat(0,G.cells.num,1); %#ok<REPMAT> % Weight Percentage
model.Mineralogy.Brucite      = repmat(0,G.cells.num,1); %#ok<REPMAT> % Weight Percentage
model.Mineralogy.Portlandite  = repmat(0,G.cells.num,1); %#ok<REPMAT> % Weight Percentage
model.Mineralogy.Pyrite       = repmat(0,G.cells.num,1); %#ok<REPMAT> % Weight Percentage
model.Mineralogy.Gypsum       = repmat(0,G.cells.num,1); %#ok<REPMAT> % Weight Percentage

%% Kinetic Data 
%Methanogenesis
model.Kinetic.mu_MET   = 1.109;%1.109;  0.3 4.1      %per day
model.Kinetic.b_MET    = 0.01 * model.Kinetic.mu_MET; %Decay Coefficent
model.Kinetic.Y_MET    = 0.03*4;       %Growth Yield
model.Kinetic.K_Dmet   = 10e-6;      %Electron Donor Half Saturation Constant
model.Kinetic.K_Amet   = 230e-6;     %Electron Acceptor Half Saturation Constant
model.Kinetic.MW_MET   = 24.6;         %Microbe Molecular Weight g/mol
model.Kinetic.m_MET    = 1e-14;        %Biomass Mass gram
model.Kinetic.N0_MET   = 1e6*1000;           %number of biomass per 1 kg of water
model.Kinetic.M0_MET   = model.Kinetic.N0_MET * model.Kinetic.m_MET * 1/model.Kinetic.MW_MET; %Initial mole of Biomass
model.Kinetic.Nmax_MET = 1e10*1000;       %Maximum biomass concnetration
model.Kinetic.Mmax_MET = model.Kinetic.Nmax_MET * model.Kinetic.m_MET * 1/model.Kinetic.MW_MET; %Maximum biomass mole per 1 kg of water
model.Kinetic.Mmin_MET = model.Kinetic.N0_MET * model.Kinetic.m_MET / model.Kinetic.MW_MET; %Minimum biomass mole per 1 kg of water
%Sulfate reduction
model.Kinetic.mu_SRB   = 1.048;%5.5;%1.048; 0.2       %per day
model.Kinetic.b_SRB    = 0.01 * model.Kinetic.mu_SRB; %Decay Coefficent
model.Kinetic.Y_SRB    = 0.08*4;       %Growth Yield
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
model.Kinetic.mu_ACE   = 0.872;%1.9;%0.872; 0.4       %per day
model.Kinetic.b_ACE    = 0.01 * model.Kinetic.mu_ACE; %Decay Coefficent
model.Kinetic.Y_ACE    = 0.07*4;       %Growth Yield
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
model.Kinetic.Y_FRB    = 0.14;       %Growth Yield
model.Kinetic.K_Dfrb   = 1e-6;       %Electron Donor Half Saturation Constant
model.Kinetic.K_Afrb   = 1e-12;     %#ok<NASGU> %Electron Acceptor Half Saturation Constant
model.Kinetic.MW_FRB   = 24.6;         %Microbe Molecular Weight g/mol
model.Kinetic.m_FRB    = 1e-14;        %Biomass Mass gram
model.Kinetic.N0_FRB   = 1e6*1000;           %numnber of biomass per 1 kg of water
model.Kinetic.M0_FRB   = model.Kinetic.N0_FRB * model.Kinetic.m_FRB * 1/model.Kinetic.MW_FRB; %Initial mole of Biomass
model.Kinetic.Nmax_FRB = 1e10*1000;       %Maximum biomass concnetration
model.Kinetic.Mmax_FRB = model.Kinetic.Nmax_FRB * model.Kinetic.m_FRB * 1/model.Kinetic.MW_FRB; %Maximum biomass mole per 1 kg of water
model.Kinetic.Mmin_FRB = model.Kinetic.N0_FRB * model.Kinetic.m_FRB / model.Kinetic.MW_FRB; %Minimum biomass mole per 1 kg of water

% stand alone flash - this section is not necessary. Here is to find Z_V
eos = EquationOfStateModel([], mixture, 'Peng-Robinson');
[L, x, y, Z_L, Z_V] = standaloneFlash(p, T, [0,1,0,0,0,0], eos);

%% Storage Scenario
totTime = 250*day;
steps   = 125;
dt = totTime/steps;
dt = repmat(dt, steps, 1);
schedule = simpleSchedule(dt);
schedule.step.control(26:100)=2; % Control Value to switch to a new schedule
schedule.step.control(101:end)=3; % Control Value to switch to a new schedule

Bg = 101325/298.15*Z_V*T/p;
rate = 0.005*pv*meter^3/day/Bg; % Surface Rate
W1 = [];
W2 = [];
W3 = [];
tmp = cell(4,1);
schedule.control=struct('W',tmp,'bc',tmp,'src',tmp);
%Injection
W1 = verticalWell(W1, G, rock, 1, 1, 1,...
                'Type', 'rate', 'Val', rate, ...
                'Name', 'Injector','comp_i',[0 1],'sign',+1,'radius',0.1);
W1(1).components = [0, 1, 0, 0, 0, 0];
%Storage
W2 = verticalWell(W2, G, rock, 1, 1, 1,...
                'Type', 'rate', 'Val', 0, ...
                'Name', 'Shut-in','comp_i',[0 1],'sign',+1,'radius',0.1);
W2(1).components = [0, 1, 0, 0, 0, 0];
%Production
W3 = verticalWell(W3, G, rock, 1, 1, 1,...
                'Type', 'grat', 'Val', -rate, ...
                'Name', 'Producer','comp_i',[0.5 0.5],'sign',-1,'radius',0.1);
W3(1).components = [0, 1, 0, 0, 0, 0];
W3.lims.bhp = p;

schedule.control(1).W=W1;
schedule.control(2).W=W2;
schedule.control(3).W=W3;
% schedule.control(4).W=W4;
s  = EOSSeparator('pressure', 1*atm, 'T', 298.15); % Set conditions for surface
sg = SeparatorGroup(s);                         % Group = single separator
sg.mode = 'moles';                              % Use mole mode 
model.FacilityModel.SeparatorGroup = sg;    % Connect to reservoir model
%% Running the simulation with a visualization tool
tic
[wellSol, states, report] = simulateScheduleAD_Modified(state0, model, schedule);
toc
