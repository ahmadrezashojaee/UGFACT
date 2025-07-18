clear;clc;close all
mrstModule add compositional ad-core ad-props mrst-gui ad-blackoil

%% Import Data from Excel
% go to Excel_Data folder and set your simulation scenario in DATA.xlsx
% This Excel workflow is only for 1D - you can expand it.
filename = 'Excel_Data\DATA.xlsx'; 

% Read the Excel data
[num, txt, raw] = xlsread(filename);

% Extract the required columns
hints = raw(:,1);        % First column (user hints)
values = raw(:,2);       % Second column (parameter values)
variable_names = raw(:,3); % Third column (MATLAB variable names)

% Loop through each row and assign values dynamically
for i = 1:length(variable_names)
    var_name = variable_names{i}; % Extract variable name
    if ischar(var_name) && ~isempty(var_name) % Ensure it's a valid name
        assignin('base', var_name, values{i}); % Save variable in MATLAB workspace
    end
end
%% Gridding
G=computeGeometry(cartGrid([NX,NY,NZ],[DX,DY,DZ])); % Making a cartesian grid
%% Rock model

rock = makeRock(G, Perm*milli*darcy, phi); % Making a rock model based on Grid, perm and poro
pv = sum(G.cells.volumes.*rock.poro); % Calculating the pore volume of each cell
%% Fluid model
f = initSimpleADIFluid('phases', 'og', 'blackoil', false); % Using O as water is not supported in seprator 
 
% Read krG.xlsx
krG_data = readmatrix('Excel_Data\krG.xlsx');  % Automatically skips header if numeric data starts after header
Sg_table = krG_data(:,1);           % First column: Sg
krG_table = krG_data(:,2);          % Second column: krG

% Read krW.xlsx
krW_data = readmatrix('Excel_Data\krW.xlsx');  % Same assumption about headers
Sw_table = krW_data(:,1);           % First column: Sw
krW_table = krW_data(:,2);          % Second column: krW

f.krG = @(sg)interpTable(Sg_table,krG_table,sg); 
f.krO = @(sw)interpTable(Sw_table,krW_table,sw); % Using O as water is not supported in seprator 


mixture = TableCompositionalMixture({'Water','Hydrogen','CarbonDioxide','Methane','HydrogenSulfide','Nitrogen'},...
                                    {'Water','H2','CO2','CH4','H2S','N2'});

arg = {G, rock, f, ... % Standard arguments
mixture,... % Compositional mixture
'water', false, 'oil', true, 'gas', true,... % Water-Gas system
'liquidPhase', 'O', 'vaporPhase', 'G'}; % Water=liquid, gas=vapor % Using O as water is not supported in seprator 
% Construct models for both formulations. Same input arguments

%% Defining the model based on the Grid, Rock, fluid, and mixture.
model = GenericOverallCompositionModel(arg{:}); % Overall mole fractions model

%model = GenericNaturalVariablesModel(arg{:}); % Natural variables
model.EOSModel.PropertyModel.volumeShift = [0, 0, 0, 0, 0, 0]; % Volume shift
gravity reset on
%% Initial state
p = Pi*barsa; T = 273.15 + Temp; s = []; z = [Z_H2O, Z_H2, Z_CO2, Z_CH4, Z_N2, Z_H2S]; % p, T, s, z

state0   = initCompositionalState(model, p, T, s, z); % Initialize state
%% Input data for biochemical and geochemical reactions
model.BioGeo = true;
model.parpool = false;   % Use parallel calculations
model.parpoolCores  = 5; % Number of cores
model.BioGeoSteps  = 5;  % Number of sub steps for bio-geochemical simulations 
model.Solution.pH    = repmat(pH,G.cells.num,1);
model.Solution.Unit  = 'mol/kgw';
%Ion concentration in the given unit
model.Solution.Ca    = repmat(Ca,G.cells.num,1);
model.Solution.Cl    = repmat(Cl,G.cells.num,1);
model.Solution.Na    = repmat(Na,G.cells.num,1);
model.Solution.K     = repmat(K,G.cells.num,1); %#ok<RPMT0>
model.Solution.S6    = repmat(S6,G.cells.num,1);
model.Solution.C4    = repmat(C4,G.cells.num,1); 
model.Solution.Mg    = repmat(Mg,G.cells.num,1);
model.Solution.Fe3   = repmat(Fe3,G.cells.num,1); %#ok<REPMAT>
model.Solution.Fe2   = repmat(Fe2,G.cells.num,1); %#ok<REPMAT>
model.Solution.Acetate   = repmat(Acetate,G.cells.num,1); %#ok<REPMAT>
model.Solution.S2    = repmat(S2,G.cells.num,1); %#ok<REPMAT>;
model.Solution.Si    = repmat(0,G.cells.num,1); %#ok<REPMAT>;

model.Mineralogy.Calcite      = repmat(Calcite,G.cells.num,1); %#ok<REPMAT> %Percent
model.Mineralogy.Dolomite     = repmat(Dolomite,G.cells.num,1); %#ok<REPMAT> % %Percent
model.Mineralogy.Quartz       = repmat(Quartz,G.cells.num,1); %#ok<REPMAT> % %Percent
model.Mineralogy.Anhydrite    = repmat(Anhydrite,G.cells.num,1); %#ok<REPMAT> %Percent
model.Mineralogy.Goethite     = repmat(Goethite,G.cells.num,1); %#ok<REPMAT> %Percent
model.Mineralogy.Brucite      = repmat(Brucite,G.cells.num,1); %#ok<REPMAT> %Percent
model.Mineralogy.Portlandite  = repmat(Portlandite,G.cells.num,1); %#ok<REPMAT> %Percent
model.Mineralogy.Pyrite       = repmat(Pyrite,G.cells.num,1); %#ok<REPMAT> %Percent
model.Mineralogy.Gypsum       = repmat(Gypsum,G.cells.num,1); %#ok<REPMAT> %Percent

%% Kinetic Data 
%Methanogenesis
model.Kinetic.mu_MET   = mu_MET;%1.109;  0.3 4.2      %per day
model.Kinetic.b_MET    = b_MET; %Decay Coefficent
model.Kinetic.Y_MET    = Y_MET;       %Growth Yield
model.Kinetic.K_Dmet   = K_Dmet;      %Electron Donor Half Saturation Constant
model.Kinetic.K_Amet   = K_Amet;     %Electron Acceptor Half Saturation Constant
model.Kinetic.MW_MET   = 24.6;         %Microbe Molecular Weight g/mol
model.Kinetic.m_MET    = 1e-14;        %Biomass Mass gram
model.Kinetic.N0_MET   = 1e6*1000;           %number of biomass per 1 kg of water
model.Kinetic.M0_MET   = model.Kinetic.N0_MET * model.Kinetic.m_MET * 1/model.Kinetic.MW_MET; %Initial mole of Biomass
model.Kinetic.Nmax_MET = 1e10*1000;       %Maximum biomass concnetration
model.Kinetic.Mmax_MET = model.Kinetic.Nmax_MET * model.Kinetic.m_MET * 1/model.Kinetic.MW_MET; %Maximum biomass mole per 1 kg of water
model.Kinetic.Mmin_MET = model.Kinetic.N0_MET * model.Kinetic.m_MET / model.Kinetic.MW_MET; %Minimum biomass mole per 1 kg of water
%Sulfate reduction
model.Kinetic.mu_SRB   = mu_SRB;%5.5;%1.048; 0.2       %per day
model.Kinetic.b_SRB    = b_SRB; %Decay Coefficent
model.Kinetic.Y_SRB    = Y_SRB;       %Growth Yield
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
model.Kinetic.mu_ACE   = mu_ACE ;%1.9;%0.872; 0.4       %per day
model.Kinetic.b_ACE    = b_ACE * model.Kinetic.mu_ACE; %Decay Coefficent
model.Kinetic.Y_ACE    = Y_ACE;       %Growth Yield
model.Kinetic.K_Dace   = K_Dace;       %Electron Donor Half Saturation Constant
model.Kinetic.K_Aace   = K_Aace;     %Electron Acceptor Half Saturation Constant
model.Kinetic.MW_ACE   = 24.6;         %Microbe Molecular Weight g/mol
model.Kinetic.m_ACE    = 1e-14;        %Biomass Mass gram
model.Kinetic.N0_ACE   = 1e6*1000;           %numnber of biomass per 1 kg of water
model.Kinetic.M0_ACE   = model.Kinetic.N0_ACE * model.Kinetic.m_ACE * 1/model.Kinetic.MW_ACE; %Initial mole of Biomass
model.Kinetic.Nmax_ACE = 1e10*1000;       %Maximum biomass concnetration
model.Kinetic.Mmax_ACE = model.Kinetic.Nmax_ACE * model.Kinetic.m_ACE * 1/model.Kinetic.MW_ACE; %Maximum biomass mole per 1 kg of water
model.Kinetic.Mmin_ACE = model.Kinetic.N0_ACE * model.Kinetic.m_ACE / model.Kinetic.MW_ACE; %Minimum biomass mole per 1 kg of water
%Iron reduction
model.Kinetic.mu_FRB   = mu_FRB;%1.5;        %per day
model.Kinetic.b_FRB    = b_FRB * model.Kinetic.mu_FRB; %Decay Coefficent
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
[L, x, y, Z_L, Z_V] = standaloneFlash(p, T, z, eos);

%% Storage Scenario
totTime = (INJ_TIME + SHUT_TIME + PROD_TIME)*day;
steps   = TOT_Steps;
dt = totTime/steps;
dt = repmat(dt, steps, 1);
schedule = simpleSchedule(dt);
schedule.step.control(INJ_Steps+1:INJ_Steps+Shut_Steps)=2; % Control Value to switch to a new schedule
schedule.step.control(INJ_Steps + Shut_Steps+1:end)=3; % Control Value to switch to a new schedule

Bg = 101325/298.15*Z_V*T/p;
rate_Injection = rate_INJ/86400; % Surface Rate
rate_Production = rate_PROD/86400; % Surface Rate
W1 = [];
W2 = [];
W3 = [];
tmp = cell(3,1);
schedule.control=struct('W',tmp,'bc',tmp,'src',tmp);
%Injection
W1 = verticalWell(W1, G, rock, 1, 1, 1,...
                'Type', 'rate', 'Val', rate_Injection, ...
                'Name', 'Injector','comp_i',[0 1],'sign',+1,'radius',0.1);
W1(1).components = [0, 1, 0, 0, 0, 0];
%Storage
W2 = verticalWell(W2, G, rock, 1, 1, 1,...
                'Type', 'rate', 'Val', 0, ...
                'Name', 'Shut-in','sign',+1,'comp_i',[0 1],'radius',0.1);
W2(1).components = [0, 1, 0, 0, 0, 0];
%Production
W3 = verticalWell(W3, G, rock, 1, 1, 1,...
                'Type', 'rate', 'Val', -rate_Production, ...
                'Name', 'Producer','comp_i',[0 1],'sign',-1,'radius',0.1);
%W3.lims.bhp = p ;
W3(1).components = [0, 1, 0, 0, 0, 0];
s  = EOSSeparator('pressure', 1*atm, 'T', 298.15); % Set conditions for surface
sg = SeparatorGroup(s);                         % Group = single separator
sg.mode = 'moles';                              % Use mole mode 
model.FacilityModel.SeparatorGroup = sg;    % Connect to reservoir model
schedule.control(1).W=W1;
schedule.control(2).W=W2;
schedule.control(3).W=W3;
%% Running the simulation with a visualization tool
tic
% fn = getPlotAfterStep(state0, model, schedule, ...
%                       'plotReservoir', true, 'view', [0, 1], ...
%                       'field', 'y:2');
[wellSol, states, report] = simulateScheduleAD_Modified(state0, model, schedule);%,'afterStepFn', fn);
toc
figure
plotToolbar(G,states,'pauseTime',0.08)