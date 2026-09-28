%clear; clc;

%% ------------------- Load simulation results ----------------
%load('LowRate.mat');      statesL = states;
%load('ModerateRate.mat'); statesM = states;
%load('HighRate.mat');     statesH = states;

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
cyclesGroup1 = [5 10 15];
cyclesGroup2 = [20 25 30];

%% ============================================================
% A) H2 in place (mol/m^3) at end of shut-in
%% ============================================================

% --- Global colour limits across all six cycles ---
allVals = [];
allCycles = [cyclesGroup1 cyclesGroup2];

for c = allCycles
    t = (c-1)*nPerCycle + idxEndShutin;

    H2L = statesL{t}.FlowProps.ComponentTotalMass{2,1} ./ MW.H2 ./ G.cells.volumes;
    H2M = statesM{t}.FlowProps.ComponentTotalMass{2,1} ./ MW.H2 ./ G.cells.volumes;
    H2H = statesH{t}.FlowProps.ComponentTotalMass{2,1} ./ MW.H2 ./ G.cells.volumes;

    allVals = [allVals; H2L; H2M; H2H];
end

climH2 = [min(allVals) max(allVals)];

%% -------- Figure 1: cycles 5, 10, 15 --------
plotH2Triplet(statesL, statesM, statesH, MW, G, X, Z, NX, NY, NZ, ...
    nPerCycle, idxEndShutin, cyclesGroup1, climH2, ...
    '');

%% -------- Figure 2: cycles 20, 25, 30 --------
plotH2Triplet(statesL, statesM, statesH, MW, G, X, Z, NX, NY, NZ, ...
    nPerCycle, idxEndShutin, cyclesGroup2, climH2, ...
    '');

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

% --- cumulative over cycles [NC x nCycles] ---
H2lossCum_L = cumsum(H2cons_L,2);
H2lossCum_M = cumsum(H2cons_M,2);
H2lossCum_H = cumsum(H2cons_H,2);

% --- Global colour limits across all six cycles ---
allLoss = [];
for c = allCycles
    lossL = H2lossCum_L(:,c) ./ G.cells.volumes;
    lossM = H2lossCum_M(:,c) ./ G.cells.volumes;
    lossH = H2lossCum_H(:,c) ./ G.cells.volumes;
    allLoss = [allLoss; lossL; lossM; lossH];
end
climLoss = [min(allLoss) max(allLoss)];

%% -------- Figure 3: cycles 5, 10, 15 --------
plotLossTriplet(H2lossCum_L, H2lossCum_M, H2lossCum_H, G, X, Z, NX, NY, NZ, ...
    cyclesGroup1, climLoss, ...
    '');

%% -------- Figure 4: cycles 20, 25, 30 --------
plotLossTriplet(H2lossCum_L, H2lossCum_M, H2lossCum_H, G, X, Z, NX, NY, NZ, ...
    cyclesGroup2, climLoss, ...
    '');

%% ============================================================
% Local functions
%% ============================================================

function plotH2Triplet(statesL, statesM, statesH, MW, G, X, Z, NX, NY, NZ, ...
                       nPerCycle, idxEndShutin, cyclesToPlot, climH2, figTitle)

    figure('Color','w');
    set(gcf,'Units','inches','Position',[1 1 18 10]);
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
        ax1 = subplot(3,3,(i-1)*3 + 1);
        surf(ax1,X,Z,H2L2D,'EdgeColor','none');
        view(ax1,2); shading(ax1,'interp');
        set(ax1,'YDir','reverse','FontSize',10,'LineWidth',0.8);
        axis(ax1,'tight');
        caxis(ax1,climH2);
        title(ax1,sprintf('Cycle %d: Low-Rate',c),'FontSize',12,'FontWeight','bold');
        xlabel(ax1,'X (m)'); ylabel(ax1,'Z (m)');
        addBigColorbar(ax1);

        % ----- Column 2: Moderate-Rate -----
        ax2 = subplot(3,3,(i-1)*3 + 2);
        surf(ax2,X,Z,H2M2D,'EdgeColor','none');
        view(ax2,2); shading(ax2,'interp');
        set(ax2,'YDir','reverse','FontSize',10,'LineWidth',0.8);
        axis(ax2,'tight');
        caxis(ax2,climH2);
        title(ax2,sprintf('Cycle %d: Moderate-Rate',c),'FontSize',12,'FontWeight','bold');
        xlabel(ax2,'X (m)'); ylabel(ax2,'Z (m)');
        addBigColorbar(ax2);

        % ----- Column 3: High-Rate -----
        ax3 = subplot(3,3,(i-1)*3 + 3);
        surf(ax3,X,Z,H2H2D,'EdgeColor','none');
        view(ax3,2); shading(ax3,'interp');
        set(ax3,'YDir','reverse','FontSize',10,'LineWidth',0.8);
        axis(ax3,'tight');
        caxis(ax3,climH2);
        title(ax3,sprintf('Cycle %d: High-Rate',c),'FontSize',12,'FontWeight','bold');
        xlabel(ax3,'X (m)'); ylabel(ax3,'Z (m)');
        addBigColorbar(ax3);
    end

    sgtitle(figTitle,'FontSize',16,'FontWeight','bold');
end

function plotLossTriplet(H2lossCum_L, H2lossCum_M, H2lossCum_H, G, X, Z, NX, NY, NZ, ...
                         cyclesToPlot, climLoss, figTitle)

    figure('Color','w');
    set(gcf,'Units','inches','Position',[1 1 18 10]);
    colormap(jet);

    for i = 1:numel(cyclesToPlot)
        c = cyclesToPlot(i);

        lossL = H2lossCum_L(:,c) ./ G.cells.volumes;
        lossM = H2lossCum_M(:,c) ./ G.cells.volumes;
        lossH = H2lossCum_H(:,c) ./ G.cells.volumes;

        lossL2D = squeeze(reshape(lossL,[NX,NY,NZ]))';
        lossM2D = squeeze(reshape(lossM,[NX,NY,NZ]))';
        lossH2D = squeeze(reshape(lossH,[NX,NY,NZ]))';

        % ----- Column 1: Low-Rate -----
        ax1 = subplot(3,3,(i-1)*3 + 1);
        surf(ax1,X,Z,lossL2D,'EdgeColor','none');
        view(ax1,2); shading(ax1,'interp');
        set(ax1,'YDir','reverse','FontSize',10,'LineWidth',0.8);
        axis(ax1,'tight');
        caxis(ax1,climLoss);
        title(ax1,sprintf('Cycle %d: Low-Rate',c),'FontSize',12,'FontWeight','bold');
        xlabel(ax1,'X (m)'); ylabel(ax1,'Z (m)');
        addBigColorbar(ax1);

        % ----- Column 2: Moderate-Rate -----
        ax2 = subplot(3,3,(i-1)*3 + 2);
        surf(ax2,X,Z,lossM2D,'EdgeColor','none');
        view(ax2,2); shading(ax2,'interp');
        set(ax2,'YDir','reverse','FontSize',10,'LineWidth',0.8);
        axis(ax2,'tight');
        caxis(ax2,climLoss);
        title(ax2,sprintf('Cycle %d: Moderate-Rate',c),'FontSize',12,'FontWeight','bold');
        xlabel(ax2,'X (m)'); ylabel(ax2,'Z (m)');
        addBigColorbar(ax2);

        % ----- Column 3: High-Rate -----
        ax3 = subplot(3,3,(i-1)*3 + 3);
        surf(ax3,X,Z,lossH2D,'EdgeColor','none');
        view(ax3,2); shading(ax3,'interp');
        set(ax3,'YDir','reverse','FontSize',10,'LineWidth',0.8);
        axis(ax3,'tight');
        caxis(ax3,climLoss);
        title(ax3,sprintf('Cycle %d: High-Rate',c),'FontSize',12,'FontWeight','bold');
        xlabel(ax3,'X (m)'); ylabel(ax3,'Z (m)');
        addBigColorbar(ax3);
    end

    sgtitle(figTitle,'FontSize',16,'FontWeight','bold');
end

function addBigColorbar(ax)
    ax.Units = 'normalized';
    axPos = ax.Position;

    % shrink axes slightly to make room
    ax.Position = [axPos(1), axPos(2), axPos(3)*0.86, axPos(4)];

    cb = colorbar(ax,'eastoutside');
    cb.FontSize = 10;
    cb.Units = 'normalized';

    axPos = ax.Position;
    cb.Position = [axPos(1)+axPos(3)+0.006, ...
                   axPos(2)+0.03*axPos(4), ...
                   0.012, ...
                   0.94*axPos(4)];
end