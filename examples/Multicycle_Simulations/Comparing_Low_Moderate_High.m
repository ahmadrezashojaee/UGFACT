%% ============================================================
%   Compare Low-, Moderate-, and High-Rate Cases
%   A) H2 in place (mol/m^3) at end of shut-in
%   B) Cumulative H2 loss per grid block (mol/m^3)
%   C) Total H2 consumption per cycle (Tonne/cycle)
%   D) Recovery and Purity vs Cycle
%% ============================================================
clear; clc;

%% ------------------- Load simulation results ----------------
load('LowRate.mat');      statesL = states;
load('ModerateRate.mat'); statesM = states;
load('HighRate.mat');     statesH = states;

load('stepDates.mat');
load('schedule.mat');

% Molecular weights [kg/mol]
MW.H2  = 2.016e-3;

%% ------------------- Load grid & geometry -------------------
fn = 'grid.GRDECL';
G  = readGRDECL(fn);
G  = processGRDECL(G);
G  = computeGeometry(G);

NX = G.cartDims(1);
NY = G.cartDims(2);
NZ = G.cartDims(3);
NC = NX*NY*NZ;

xc = reshape(G.cells.centroids(:,1), G.cartDims);
zc = reshape(G.cells.centroids(:,3), G.cartDims);

X = squeeze(xc(:,1,:))';    % [NZ, NX]
Z = squeeze(zc(:,1,:))';    % [NZ, NX]

%% ---------------------- Safety check ------------------------
if size(statesL{1}.y,1) ~= NC
    error('Grid / state size mismatch with LowRate.');
end
if size(statesM{1}.y,1) ~= NC || size(statesH{1}.y,1) ~= NC
    error('Grid / state size mismatch in ModerateRate/HighRate.');
end

%% -------------------- Cycle indexing ------------------------
nPerCycle    = 75;          % 31 inj + 13 shut-in + 31 prod
idxEndShutin = 44;          % end of shut-in within each cycle
cyclesToPlot = [5 10 15 20 25 30];

%% ============================================================
% A) H2 in place (mol/m^3) at end of shut-in
%% ============================================================

allVals = [];

for c = cyclesToPlot
    t = (c-1)*nPerCycle + idxEndShutin;

    H2L = statesL{t}.FlowProps.ComponentTotalMass{2,1} ./ MW.H2 ./ G.cells.volumes;
    H2M = statesM{t}.FlowProps.ComponentTotalMass{2,1} ./ MW.H2 ./ G.cells.volumes;
    H2H = statesH{t}.FlowProps.ComponentTotalMass{2,1} ./ MW.H2 ./ G.cells.volumes;

    allVals = [allVals; H2L; H2M; H2H];
end

climH2 = [min(allVals) max(allVals)];

figure('Color','w');
set(gcf,'Units','normalized','Position',[0 0 1 1]);
tiledlayout(6,3,'TileSpacing','tight','Padding','compact');
colormap(jet);

for i = 1:numel(cyclesToPlot)
    c = cyclesToPlot(i);
    t = (c-1)*nPerCycle + idxEndShutin;

    H2L = statesL{t}.FlowProps.ComponentTotalMass{2,1} ./ MW.H2 ./ G.cells.volumes;
    H2M = statesM{t}.FlowProps.ComponentTotalMass{2,1} ./ MW.H2 ./ G.cells.volumes;
    H2H = statesH{t}.FlowProps.ComponentTotalMass{2,1} ./ MW.H2 ./ G.cells.volumes;

    H2L2D = squeeze(reshape(H2L,[NX,NY,NZ]))';
    H2M2D = squeeze(reshape(H2M,[NX,NY,NZ]))';
    H2H2D = squeeze(reshape(H2H,[NX,NY,NZ]))';

    % ----- Column 1: Low-Rate -----
    nexttile;
    surf(X,Z,H2L2D,'EdgeColor','none');
    view(2); shading interp; set(gca,'YDir','reverse'); axis tight;
    caxis(climH2);
    title(sprintf('Cycle %d: Low-Rate',c),'FontSize',10);
    xlabel('X (m)'); ylabel('Z (m)');
    colorbar;

    % ----- Column 2: Moderate-Rate -----
    nexttile;
    surf(X,Z,H2M2D,'EdgeColor','none');
    view(2); shading interp; set(gca,'YDir','reverse'); axis tight;
    caxis(climH2);
    title(sprintf('Cycle %d: Moderate-Rate',c),'FontSize',10);
    xlabel('X (m)'); ylabel('Z (m)');
    colorbar;

    % ----- Column 3: High-Rate -----
    nexttile;
    surf(X,Z,H2H2D,'EdgeColor','none');
    view(2); shading interp; set(gca,'YDir','reverse'); axis tight;
    caxis(climH2);
    title(sprintf('Cycle %d: High-Rate',c),'FontSize',10);
    xlabel('X (m)'); ylabel('Z (m)');
    colorbar;
end

sgtitle('H_2 in place (mol/m^3) at the end of shut-in for selected cycles', ...
        'FontSize',14,'FontWeight','bold');

%% ============================================================
% B) Cumulative H2 loss per grid block (mol/m^3)
%    Computed from MET + SRB + ACE rates
%% ============================================================

H2cons_L = computeH2consCycles(statesL, stepDates, nPerCycle);
H2cons_M = computeH2consCycles(statesM, stepDates, nPerCycle);
H2cons_H = computeH2consCycles(statesH, stepDates, nPerCycle);

[NCcheck, nCycles] = size(H2cons_L);
if NCcheck ~= NC
    error('H2cons_L has wrong number of cells.');
end

% Cumulative over cycles [NC x nCycles]
H2lossCum_L = cumsum(H2cons_L,2);
H2lossCum_M = cumsum(H2cons_M,2);
H2lossCum_H = cumsum(H2cons_H,2);

% Global colour limits across selected cycles
allLoss = [];
for c = cyclesToPlot
    lossL = H2lossCum_L(:,c) ./ G.cells.volumes;
    lossM = H2lossCum_M(:,c) ./ G.cells.volumes;
    lossH = H2lossCum_H(:,c) ./ G.cells.volumes;
    allLoss = [allLoss; lossL; lossM; lossH];
end
climLoss = [min(allLoss) max(allLoss)];

figure('Color','w');
set(gcf,'Units','normalized','Position',[0 0 1 1]);
tiledlayout(6,3,'TileSpacing','tight','Padding','compact');
colormap(jet);

for i = 1:numel(cyclesToPlot)
    c = cyclesToPlot(i);

    lossL = H2lossCum_L(:,c) ./ G.cells.volumes;   % mol/m^3
    lossM = H2lossCum_M(:,c) ./ G.cells.volumes;
    lossH = H2lossCum_H(:,c) ./ G.cells.volumes;

    lossL2D = squeeze(reshape(lossL,[NX,NY,NZ]))';
    lossM2D = squeeze(reshape(lossM,[NX,NY,NZ]))';
    lossH2D = squeeze(reshape(lossH,[NX,NY,NZ]))';

    % ----- Low-Rate -----
    nexttile;
    surf(X,Z,lossL2D,'EdgeColor','none');
    view(2); shading interp; set(gca,'YDir','reverse'); axis tight;
    caxis(climLoss);
    title(sprintf('Cycle %d: Low-Rate',c),'FontSize',10);
    xlabel('X (m)'); ylabel('Z (m)');
    colorbar;

    % ----- Moderate-Rate -----
    nexttile;
    surf(X,Z,lossM2D,'EdgeColor','none');
    view(2); shading interp; set(gca,'YDir','reverse'); axis tight;
    caxis(climLoss);
    title(sprintf('Cycle %d: Moderate-Rate',c),'FontSize',10);
    xlabel('X (m)'); ylabel('Z (m)');
    colorbar;

    % ----- High-Rate -----
    nexttile;
    surf(X,Z,lossH2D,'EdgeColor','none');
    view(2); shading interp; set(gca,'YDir','reverse'); axis tight;
    caxis(climLoss);
    title(sprintf('Cycle %d: High-Rate',c),'FontSize',10);
    xlabel('X (m)'); ylabel('Z (m)');
    colorbar;
end

sgtitle('Cumulative H_2 loss per grid block (mol/m^3) up to each cycle', ...
        'FontSize',14,'FontWeight','bold');

%% ============================================================
% C) Total and cumulative H2 consumed per cycle (Tonne/cycle)
%% ============================================================
% Define colours
colLow  = [0 0.4470 0.7410];       % blue
colMod  = [0.8500 0.3250 0.0980];  % orange
colHigh = [0 0 0];                 % black
cycles = 1:nCycles;

H2tot_L = sum(H2cons_L,1);   % [1 x nCycles] mol/cycle
H2tot_M = sum(H2cons_M,1);
H2tot_H = sum(H2cons_H,1);

H2cum_L = cumsum(H2tot_L);   % cumulative mol
H2cum_M = cumsum(H2tot_M);
H2cum_H = cumsum(H2tot_H);

figure('Color','w'); clf;

% -------- Left panel: per-cycle consumption --------
subplot(1,2,1);
plot(cycles, H2tot_L.*MW.H2./1e3, '-o','LineWidth',1.5,'MarkerSize', 6, 'MarkerFaceColor', colLow, 'Color', colLow); hold on;
plot(cycles, H2tot_M.*MW.H2./1e3, '-s','LineWidth',1.5,'MarkerSize', 6, 'MarkerFaceColor', colMod, 'Color', colMod);
plot(cycles, H2tot_H.*MW.H2./1e3, '-^','LineWidth',1.5,'MarkerSize', 6, 'MarkerFaceColor', colHigh, 'Color', colHigh);
grid on;
xlabel('Cycle number');
ylabel('H_2 consumed per cycle [Tonne]');
legend({'Low-Rate','Moderate-Rate','High-Rate'},'Location','best');
title('H_2 consumption per cycle');

% -------- Right panel: cumulative consumption --------
subplot(1,2,2);
plot(cycles, H2cum_L.*MW.H2./1e3, '-o','LineWidth',1.5,'MarkerSize', 6, 'MarkerFaceColor', colLow, 'Color', colLow); hold on;
plot(cycles, H2cum_M.*MW.H2./1e3, '-s','LineWidth',1.5,'MarkerSize', 6, 'MarkerFaceColor', colMod, 'Color', colMod);
plot(cycles, H2cum_H.*MW.H2./1e3, '-^','LineWidth',1.5,'MarkerSize', 6, 'MarkerFaceColor', colHigh, 'Color', colHigh);
grid on;
xlabel('Cycle number');
ylabel('Cumulative H_2 consumed [Tonne]');
legend({'Low-Rate','Moderate-Rate','High-Rate'},'Location','best');
title('Cumulative H_2 consumption over cycles');

%% ============================================================
% D) Recovery and Purity vs Cycle for different rate cases
%   - LowRate (blue), ModerateRate (orange), HighRate (black)
%% ============================================================
load('HighRate_wellSol.mat')
wellSolH = wellSol;
load('ModerateRate_wellSol.mat')
wellSolM = wellSol;
load('LowRate_wellSol.mat')
wellSolL = wellSol;

% Analyze cycles for each rate case
[purityL, recoveryL, timeL, H2SL, avgH2SL, bhpL, injH2L] = ...
    analyzeCycles_Well(wellSolL, schedule, stepDates);

[purityM, recoveryM, timeM, H2SM, avgH2SM, bhpM, injH2M] = ...
    analyzeCycles_Well(wellSolM, schedule, stepDates);

[purityH, recoveryH, timeH, H2SH, avgH2SH, bhpH, injH2H] = ...
    analyzeCycles_Well(wellSolH, schedule, stepDates);

% Number of cycles (assumed same for all)
nCycles = numel(recoveryL);
cycles  = 1:nCycles;

figure('Color','w'); clf;

%% 1) Final Recovery vs cycle
finalRec_L = 100 * cellfun(@(x) x(end), recoveryL);
finalRec_M = 100 * cellfun(@(x) x(end), recoveryM);
finalRec_H = 100 * cellfun(@(x) x(end), recoveryH);

subplot(1,2,1);
hold on;
plot(cycles, finalRec_L, '-o', 'LineWidth', 1.4, ...
    'MarkerSize', 6, 'MarkerFaceColor', colLow, 'Color', colLow);
plot(cycles, finalRec_M, '-s', 'LineWidth', 1.4, ...
    'MarkerSize', 6, 'MarkerFaceColor', colMod, 'Color', colMod);
plot(cycles, finalRec_H, '-^', 'LineWidth', 1.4, ...
    'MarkerSize', 6, 'MarkerFaceColor', colHigh, 'Color', colHigh);

grid on;
xlabel('Cycle');
ylabel('Final Recovery [%]');
title('Final H_2 Recovery per Cycle');
legend({'Low Rate','Moderate Rate','High Rate'}, 'Location', 'best');
hold off;

%% 2) Average Purity vs cycle
avgPur_L = 100 * cellfun(@mean, purityL);
avgPur_M = 100 * cellfun(@mean, purityM);
avgPur_H = 100 * cellfun(@mean, purityH);

subplot(1,2,2);
hold on;
plot(cycles, avgPur_L, '-o', 'LineWidth', 1.4, ...
    'MarkerSize', 6, 'MarkerFaceColor', colLow, 'Color', colLow);
plot(cycles, avgPur_M, '-s', 'LineWidth', 1.4, ...
    'MarkerSize', 6, 'MarkerFaceColor', colMod, 'Color', colMod);
plot(cycles, avgPur_H, '-^', 'LineWidth', 1.4, ...
    'MarkerSize', 6, 'MarkerFaceColor', colHigh, 'Color', colHigh);

grid on;
xlabel('Cycle');
ylabel('Average Purity [%]');
title('Average H_2 Purity per Cycle');
legend({'Low Rate','Moderate Rate','High Rate'}, 'Location', 'best');
hold off;
