
%% --- Load simulation results ---
load("HighRate.mat");   % with reaction
states1 = states;
load("NoReaction.mat");     % no reaction
states2 = states;

load('stepDates.mat');
load('schedule.mat');

MW.H2  = 2.016e-3; 
MW.CH4 = 16.04e-3; 
MW.CO2 = 44.01e-3; 
MW.H2S = 34.08e-3;
%% --- Load grid and build geometry ---
fn = 'grid.GRDECL';
G  = readGRDECL(fn);
G  = processGRDECL(G);
G  = computeGeometry(G);

NX = G.cartDims(1);
NY = G.cartDims(2);
NZ = G.cartDims(3);

xc = reshape(G.cells.centroids(:,1), G.cartDims);
zc = reshape(G.cells.centroids(:,3), G.cartDims);

X = squeeze(xc(:,1,:))';    % [NZ, NX]
Z = squeeze(zc(:,1,:))';    % [NZ, NX]

%% --- Safety ---
NC_expected = NX*NY*NZ;
NC_actual   = size(states1{1}.y,1);
if NC_expected ~= NC_actual
    error('Grid mismatch.');
end

%% --- Time indexing ---
nPerCycle    = 75;
idxEndShutin = 44;   % 31 + 13 shut-in
cyclesToPlot = [5 10 15 20 25 30];

%% ---------------------------------------------------------
% 1) Precompute global color limits for each column
%% ---------------------------------------------------------
all_with = [];
all_nore = [];
all_diff = [];

for c = cyclesToPlot
    t = (c - 1)*nPerCycle + idxEndShutin;

    y_with = states1{t}.y(:,2);
    y_nore = states2{t}.y(:,2);

    Yw = squeeze(reshape(y_with, [NX, NY, NZ]))';
    Yn = squeeze(reshape(y_nore, [NX, NY, NZ]))';
    Yd = Yn - Yw;

    all_with = [all_with; Yw(:)];
    all_nore = [all_nore; Yn(:)];
    all_diff = [all_diff; Yd(:)];
end

clim_with = [min(all_with) max(all_with)];
clim_nore = [min(all_nore) max(all_nore)];
clim_diff = [min(all_diff) max(all_diff)];

%% ---------------------------------------------------------
%% ---------------------------------------------------------
% 2) Publication-ready figure with tight spacing
%% ---------------------------------------------------------

figure('Color','w');
set(gcf, 'Units','normalized', 'Position',[0 0 1 1]);   % maximise screen space

tiledlayout(6,3, 'TileSpacing','tight', 'Padding','compact');
colormap(jet);

for ic = 1:numel(cyclesToPlot)
    c = cyclesToPlot(ic);
    t = (c - 1)*nPerCycle + idxEndShutin;

    y_with = states1{t}.y(:,2);
    y_nore = states2{t}.y(:,2);

    Y_with2D = squeeze(reshape(y_with, [NX, NY, NZ]))';
    Y_nore2D = squeeze(reshape(y_nore, [NX, NY, NZ]))';
    Y_diff2D = Y_nore2D - Y_with2D;

    % ---------- Column 1: With Reaction ----------
    nexttile;
    surf(X, Z, Y_with2D, 'EdgeColor','none');
    view(2); shading interp; set(gca,'YDir','reverse'); axis tight;
    caxis(clim_with);
    title(sprintf('Cycle %d: Low-Rate', c), 'FontSize',10);
    xlabel('X (m)'); ylabel('Z (m)');
    hcb = colorbar; hcb.Location = 'eastoutside';

    % ---------- Column 2: No Reaction ----------
    nexttile;
    surf(X, Z, Y_nore2D, 'EdgeColor','none');
    view(2); shading interp; set(gca,'YDir','reverse'); axis tight;
    caxis(clim_nore);
    title(sprintf('Cycle %d: No-Reaction', c), 'FontSize',10);
    xlabel('X (m)'); ylabel('Z (m)');
    hcb = colorbar; hcb.Location = 'eastoutside';

    % ---------- Column 3: Difference ----------
    nexttile;
    surf(X, Z, Y_diff2D, 'EdgeColor','none');
    view(2); shading interp; set(gca,'YDir','reverse'); axis tight;
    caxis(clim_diff);
    title(sprintf('Cycle %d: Difference', c), 'FontSize',10);
    xlabel('X (m)'); ylabel('Z (m)');
    hcb = colorbar; hcb.Location = 'eastoutside';
end

sgtitle('Mole fraction of H2 in the gas phase at the end of shut-in for selected cycles', ...
        'FontSize',14, 'FontWeight','bold');


%% ---------------------------------------------------------
% 3) H2 Moles
%% ---------------------------------------------------------
all_with = [];
all_nore = [];
all_diff = [];

for c = cyclesToPlot
    t = (c - 1)*nPerCycle + idxEndShutin;

    y_with = states1{t}.FlowProps.ComponentTotalMass{2,1}./MW.H2./1000./G.cells.volumes; %Kg.mol/rm^3
    y_nore = states2{t}.FlowProps.ComponentTotalMass{2,1}./MW.H2./1000./G.cells.volumes; %Kg.mol/rm^3

    Yw = squeeze(reshape(y_with, [NX, NY, NZ]))';
    Yn = squeeze(reshape(y_nore, [NX, NY, NZ]))';
    Yd = Yn - Yw;

    all_with = [all_with; Yw(:)];
    all_nore = [all_nore; Yn(:)];
    all_diff = [all_diff; Yd(:)];
end

clim_with = [min(all_with) max(all_with)];
clim_nore = [min(all_nore) max(all_nore)];
clim_diff = [min(all_diff) max(all_diff)];

%% ---------------------------------------------------------
%% ---------------------------------------------------------
% 4) Publication-ready figure with tight spacing
%% ---------------------------------------------------------

figure('Color','w');
set(gcf, 'Units','normalized', 'Position',[0 0 1 1]);   % maximise screen space

tiledlayout(6,3, 'TileSpacing','tight', 'Padding','compact');
colormap(jet);

for ic = 1:numel(cyclesToPlot)
    c = cyclesToPlot(ic);
    t = (c - 1)*nPerCycle + idxEndShutin;

    y_with = states1{t}.FlowProps.ComponentTotalMass{2,1}./MW.H2./1000./G.cells.volumes; %Kg.mol/rm^3
    y_nore = states2{t}.FlowProps.ComponentTotalMass{2,1}./MW.H2./1000./G.cells.volumes; %Kg.mol/rm^3

    Y_with2D = squeeze(reshape(y_with, [NX, NY, NZ]))';
    Y_nore2D = squeeze(reshape(y_nore, [NX, NY, NZ]))';
    Y_diff2D = Y_nore2D - Y_with2D;

    % ---------- Column 1: With Reaction ----------
    nexttile;
    surf(X, Z, Y_with2D, 'EdgeColor','none');
    view(2); shading interp; set(gca,'YDir','reverse'); axis tight;
    caxis(clim_nore);
    title(sprintf('Cycle %d: Low-Rate', c), 'FontSize',10);
    xlabel('X (m)'); ylabel('Z (m)');
    hcb = colorbar; hcb.Location = 'eastoutside';

    % ---------- Column 2: No Reaction ----------
    nexttile;
    surf(X, Z, Y_nore2D, 'EdgeColor','none');
    view(2); shading interp; set(gca,'YDir','reverse'); axis tight;
    caxis(clim_nore);
    title(sprintf('Cycle %d: No-Reaction', c), 'FontSize',10);
    xlabel('X (m)'); ylabel('Z (m)');
    hcb = colorbar; hcb.Location = 'eastoutside';

    % ---------- Column 3: Difference ----------
    nexttile;
    surf(X, Z, Y_diff2D, 'EdgeColor','none');
    view(2); shading interp; set(gca,'YDir','reverse'); axis tight;
    caxis(clim_diff);
    title(sprintf('Cycle %d: Difference', c), 'FontSize',10);
    xlabel('X (m)'); ylabel('Z (m)');
    hcb = colorbar; hcb.Location = 'eastoutside';
end

sgtitle('The amount of H2 (kg.mol) per each grid block volume (m^3) at the end of shut-in for selected cycles', ...
        'FontSize',14, 'FontWeight','bold');


%% 5) Dolomite
%% ---------------------------------------------------------
%% Dolomite and Calcite content per cycle + spatial maps
figure; clf;
set(gcf,'Color','w');
colormap(turbo);

% --- Find end index of each cycle ---
ctrl = schedule.step.control;
cycleEndIdx = find(ctrl(1:end-1)==3 & ctrl(2:end)==1);
cycleEndIdx = [cycleEndIdx; numel(ctrl)];
nCycles = numel(cycleEndIdx);

%% ----------------------------------------------------------
% 1) DOLOMITE — per-cycle trends
%% ----------------------------------------------------------
dolTop  = zeros(nCycles,1);
dolBot  = zeros(nCycles,1);

for c = 1:nCycles
    tEnd = cycleEndIdx(c);
    dol  = states1{tEnd}.Mineralogy.Dolomite;    % length 500

    dolTop(c) = sum(dol(1:250)   ./ G.cells.volumes(1:250));
    dolBot(c) = sum(dol(251:500) ./ G.cells.volumes(251:500));
end

subplot(2,2,1);
yyaxis left
plot(1:nCycles, dolTop, '-o', 'LineWidth',1.5);
ylabel('Dolomite in top layers (mol/m^3)');

yyaxis right
plot(1:nCycles, dolBot, '-s','LineWidth',1.5);
ylabel('Dolomite in bottom layers (mol/m^3)');

xlabel('Cycle number');
title('Dolomite content per cycle');
grid on;
legend({'Top Layers (Dolomite-Rich)', 'Bottom Layers'}, 'Location','best');


%% ----------------------------------------------------------
% 2) DOLOMITE — final 2D distribution
%% ----------------------------------------------------------
subplot(2,2,2);

tEnd = cycleEndIdx(end);     % final step index
dolEnd = states1{tEnd}.Mineralogy.Dolomite ./ G.cells.volumes;

Dol3D = reshape(dolEnd, [NX, NY, NZ]);
Dol2D = squeeze(Dol3D(:,1,:))';

surf(X, Z, Dol2D, 'EdgeColor','none');
view(2); shading interp; axis tight;
set(gca,'YDir','reverse');
xlabel('X (m)'); ylabel('Z (m)');
title('Dolomite distribution at end of simulation (mol/m^3)');
colorbar;


%% ----------------------------------------------------------
% 3) CALCITE — per-cycle trends
%% ----------------------------------------------------------
calcTop = zeros(nCycles,1);
calcBot = zeros(nCycles,1);

for c = 1:nCycles
    tEnd = cycleEndIdx(c);
    calc = states1{tEnd}.Mineralogy.Calcite;    % length 500

    calcTop(c) = sum(calc(1:250)   ./ G.cells.volumes(1:250));
    calcBot(c) = sum(calc(251:500) ./ G.cells.volumes(251:500));
end

subplot(2,2,3);
yyaxis left
plot(1:nCycles, calcTop, '-o','LineWidth',1.5);
ylabel('Calcite in top layers (mol/m^3)');

yyaxis right
plot(1:nCycles, calcBot, '-s','LineWidth',1.5);
ylabel('Calcite in bottom layers (mol/m^3)');

xlabel('Cycle number');
title('Calcite content per cycle');
grid on;
legend({'Top Layers (Dolomite-Rich)', 'Bottom Layers'}, 'Location','best');


%% ----------------------------------------------------------
% 4) CALCITE — final 2D distribution
%% ----------------------------------------------------------
subplot(2,2,4);

calcEnd = states1{tEnd}.Mineralogy.Calcite ./ G.cells.volumes;

Calc3D = reshape(calcEnd, [NX, NY, NZ]);
Calc2D = squeeze(Calc3D(:,1,:))';

surf(X, Z, Calc2D, 'EdgeColor','none');
view(2); shading interp; axis tight;
set(gca,'YDir','reverse');
xlabel('X (m)'); ylabel('Z (m)');
title('Calcite distribution at end of simulation (mol/m^3)');
colorbar;

%% ============================================================
%   CYCLE-AVERAGED MICROBIAL RATES (per kg water AND total)
%   This script uses states1 (reactive case), schedule, stepDates
%% ============================================================

statesR = states1;
nSteps  = numel(statesR);
nPerCycle = 75;                 % 31 inj + 13 shut-in + 31 prod
nCycles = nSteps / nPerCycle;   % should be 30

% Preallocate
avgMET_mass   = zeros(nCycles,1);
avgSRB_mass   = zeros(nCycles,1);
avgACE_mass   = zeros(nCycles,1);
avgTOTAL_mass = zeros(nCycles,1);

avgMET   = zeros(nCycles,1);    % mol/day
avgSRB   = zeros(nCycles,1);
avgACE   = zeros(nCycles,1);
avgTOTAL = zeros(nCycles,1);

for c = 1:nCycles
    idxRange = (c-1)*nPerCycle + (1:nPerCycle);  % timesteps in this cycle

    % Mass-normalised accumulators
    numMET_mass = 0; den_mass = 0;
    numSRB_mass = 0;
    numACE_mass = 0;

    % Total H2 consumption accumulators (mol/day)
    MET_tot = 0; SRB_tot = 0; ACE_tot = 0;
    Tcyc = 0;   % cycle duration in days

    % Loop over timesteps of the cycle
    for k = 1:numel(idxRange)
        t = idxRange(k);
        
        % time increment (days)
        if t == 1
            dt = 0;
        else
            dt = days(stepDates(t) - stepDates(t-1));
        end
        if dt == 0
            continue;
        end
        Tcyc = Tcyc + dt;

        st = statesR{t,1};
        if ~isfield(st, 'Solution')
            continue;
        end

        % Rates [mmol/day/kgw]
        r_MET = st.Solution.MET_Rate;
        r_SRB = st.Solution.SRB_Rate;
        r_ACE = st.Solution.ACE_Rate;

        % Water mass per cell [kg]
        m_water = st.FlowProps.ComponentTotalMass{1,1};

        % Active cells = cells where rate is nonzero for that process
        active_MET = (r_MET ~= 0);
        active_SRB = (r_SRB ~= 0);
        active_ACE = (r_ACE ~= 0);

        %% ----------- Mass-normalised rates -----------
        if any(active_MET)
            numMET_mass = numMET_mass + sum(r_MET(active_MET).*m_water(active_MET)) * dt;
            den_mass    = den_mass    + sum(m_water(active_MET)) * dt;
        end
        if any(active_SRB)
            numSRB_mass = numSRB_mass + sum(r_SRB(active_SRB).*m_water(active_SRB)) * dt;
        end
        if any(active_ACE)
            numACE_mass = numACE_mass + sum(r_ACE(active_ACE).*m_water(active_ACE)) * dt;
        end

        %% ----------- Total consumption rates (mol/day) -----------
        % Convert mmol/day/kgw * kgw = mmol/day → mol/day
        MET_tot = MET_tot + sum(r_MET(active_MET).*m_water(active_MET)) * dt * 1e-3;
        SRB_tot = SRB_tot + sum(r_SRB(active_SRB).*m_water(active_SRB)) * dt * 1e-3;
        ACE_tot = ACE_tot + sum(r_ACE(active_ACE).*m_water(active_ACE)) * dt * 1e-3;
    end

    % Final averages for this cycle
    avgMET_mass(c)   = numMET_mass / max(den_mass, eps);
    avgSRB_mass(c)   = numSRB_mass / max(den_mass, eps);
    avgACE_mass(c)   = numACE_mass / max(den_mass, eps);
    avgTOTAL_mass(c) = (numMET_mass + numSRB_mass + numACE_mass) / max(den_mass, eps);

    avgMET(c)   = MET_tot   / max(Tcyc, eps);
    avgSRB(c)   = SRB_tot   / max(Tcyc, eps);
    avgACE(c)   = ACE_tot   / max(Tcyc, eps);
    avgTOTAL(c) = (MET_tot + SRB_tot + ACE_tot) / max(Tcyc, eps);
end

%% ============================================================
%   PLOTTING
%% ============================================================

figure; clf; set(gcf,'Color','w');

%% ---------------- LEFT: per kg water ----------------
subplot(1,2,1)
plot(1:nCycles, avgMET_mass,   '-o', 'LineWidth',1.5); hold on;
plot(1:nCycles, avgSRB_mass,   '-s', 'LineWidth',1.5);
plot(1:nCycles, avgACE_mass,   '-^', 'LineWidth',1.5);
plot(1:nCycles, avgTOTAL_mass, '-d', 'LineWidth',1.5);

xlabel('Cycle number')
ylabel('Mass-normalised rate [mmol/day/kg_w]')
title('Microbial reaction intensity (per kg of water)')
grid on
legend({'MET','SRB','ACE','Total'}, 'Location','best')


%% ---------------- RIGHT: total H2 consumption ----------------
subplot(1,2,2)
plot(1:nCycles, avgMET,   '-o', 'LineWidth',1.5); hold on;
plot(1:nCycles, avgSRB,   '-s', 'LineWidth',1.5);
plot(1:nCycles, avgACE,   '-^', 'LineWidth',1.5);
plot(1:nCycles, avgTOTAL, '-d', 'LineWidth',1.5);

xlabel('Cycle number')
ylabel('Total H_2 consumption rate [mol/day]')
title('Total rate (mol/day) in each cycle')
grid on
legend({'MET','SRB','ACE','Total'}, 'Location','best')


%% === pH distribution at end of selected cycles (dry cells masked) ===

statesR = states1;   % or states1 if that is your reactive case

% --- Cycle end indices from schedule (end of production) ---
ctrl = schedule.step.control;
cycleEndIdx = find(ctrl(1:end-1) == 3 & ctrl(2:end) == 1);   % 3 -> 1
cycleEndIdx = [cycleEndIdx; numel(ctrl)];                     % last cycle
nCycles = numel(cycleEndIdx);

cyclesToPlot = [1 2 3 4 5 6];

% --- Precompute global pH limits over selected cycles (ignore dry NaNs) ---
allpH = [];

for c = cyclesToPlot
    tEnd = cycleEndIdx(c);

    pHvec = statesR{tEnd,1}.Solution.pH;    % vector [NC x 1]
    Sw    = statesR{tEnd,1}.s(:,1);         % water saturation

    % Mask dry cells (no pH meaning)
    pHvec(Sw < 1e-5) = NaN;

    pH3D = reshape(pHvec, [NX, NY, NZ]);
    pH2D = squeeze(pH3D(:,1,:))';          % [NZ, NX]

    allpH = [allpH; pH2D(~isnan(pH2D))];
end

clim_pH = [min(allpH) max(allpH)];

%% --- Plot pH for selected cycles ---
figure; clf;
set(gcf,'Color','w');
colormap(jet);

for i = 1:numel(cyclesToPlot)
    c    = cyclesToPlot(i);
    tEnd = cycleEndIdx(c);

    % Extract pH and Sw
    pHvec = statesR{tEnd,1}.Solution.pH;
    Sw    = statesR{tEnd,1}.s(:,1);

    % Mask dry-out cells
    pHvec(Sw < 1e-4) = NaN;

    % Reshape to 2D [NZ, NX]
    pH3D = reshape(pHvec, [NX, NY, NZ]);
    pH2D = squeeze(pH3D(:,1,:))';          % [NZ, NX]

    subplot(2,3,i);
    surf(X, Z, pH2D, 'EdgeColor','none');
    view(2); shading interp;
    set(gca,'YDir','reverse');
    axis tight;
    caxis(clim_pH);

    xlabel('X (m)');
    ylabel('Z (m)');
    title(sprintf('Cycle %d: pH', c));
    colorbar;
end

%% === MET-RATE distribution at end of selected cycles (dry cells masked) ===

statesR = states1;   % or states1 if that is your reactive case

% --- Cycle end indices from schedule (end of production) ---
ctrl = schedule.step.control;
cycleEndIdx = find(ctrl(1:end-1) == 3 & ctrl(2:end) == 1);   % 3 -> 1
cycleEndIdx = [cycleEndIdx; numel(ctrl)];                     % last cycle
nCycles = numel(cycleEndIdx);

cyclesToPlot = [1 2 3 4 5 30];

% --- Precompute global pH limits over selected cycles (ignore dry NaNs) ---
allpH = [];

for c = cyclesToPlot
    tEnd = cycleEndIdx(c);

    pHvec = statesR{tEnd,1}.Solution.MET_Rate;    % vector [NC x 1]
    Sw    = statesR{tEnd,1}.s(:,1);         % water saturation

    % Mask dry cells (no pH meaning)
    pHvec(Sw < 1e-5) = NaN;

    pH3D = reshape(pHvec, [NX, NY, NZ]);
    pH2D = squeeze(pH3D(:,1,:))';          % [NZ, NX]

    allpH = [allpH; pH2D(~isnan(pH2D))];
end

clim_pH = [min(allpH) max(allpH)];

%% --- Plot pH for selected cycles ---
figure; clf;
set(gcf,'Color','w');
colormap(jet);

for i = 1:numel(cyclesToPlot)
    c    = cyclesToPlot(i);
    tEnd = cycleEndIdx(c);

    % Extract pH and Sw
    pHvec = statesR{tEnd,1}.Solution.MET_Rate;
    Sw    = statesR{tEnd,1}.s(:,1);

    % Mask dry-out cells
    pHvec(Sw < 1e-4) = NaN;

    % Reshape to 2D [NZ, NX]
    pH3D = reshape(pHvec, [NX, NY, NZ]);
    pH2D = squeeze(pH3D(:,1,:))';          % [NZ, NX]

    subplot(2,3,i);
    surf(X, Z, pH2D, 'EdgeColor','none');
    view(2); shading interp;
    set(gca,'YDir','reverse');
    axis tight;
    caxis(clim_pH);

    xlabel('X (m)');
    ylabel('Z (m)');
    title(sprintf('Cycle %d: MET Rate (mmol/day/kg_w)', c));
    colorbar;
end

