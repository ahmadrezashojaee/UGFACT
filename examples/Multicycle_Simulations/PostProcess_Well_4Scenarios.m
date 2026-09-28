%% Load data
load('ModerateRate_wellSol.mat');  wellSol_mod = wellSol;
load('NoReaction_WS.mat');         wellSol_nr  = wellSol;
load('HighRate_wellSol.mat');      wellSol_high = wellSol;
load('LowRate_wellSol.mat');       wellSol_low  = wellSol;

load('stepDates.mat');
load('schedule.mat');

% --- Analyse all datasets ---
[pur_nr, rec_nr, time_nr, H2S_nr, avgH2S_nr, bhp_nr, inj_nr] = analyzeCycles_Well(wellSol_nr,  schedule, stepDates);
[pur_low,rec_low,time_low,H2S_low,avgH2S_low,bhp_low,inj_low] = analyzeCycles_Well(wellSol_low, schedule, stepDates);
[pur_mod,rec_mod,time_mod,H2S_mod,avgH2S_mod,bhp_mod,inj_mod] = analyzeCycles_Well(wellSol_mod, schedule, stepDates);
[pur_high,rec_high,time_high,H2S_high,avgH2S_high,bhp_high,inj_high] = analyzeCycles_Well(wellSol_high, schedule, stepDates);

nCycles = numel(pur_nr);


%% ============================================================
% 1) Recovery vs time – all cycles (for all scenarios)
%    (Only plot moderate here like your original code)
% ============================================================
figure; hold on;
for c = 1:nCycles
    plot(time_mod{c}, 100*rec_mod{c}, 'LineWidth', 1.2, ...
        'DisplayName', sprintf('Cycle %d', c));
end
xlabel('Date'); ylabel('H_2 Recovery [%]');
title('H_2 Recovery vs Time (Moderate Rate)');
legend show; grid on;



%% ============================================================
% 2) Final Recovery vs cycle (lines + grouped differences)
%    Order: No → Low → Moderate → High
% ============================================================
final_nr   = 100*cellfun(@(x) x(end), rec_nr);
final_low  = 100*cellfun(@(x) x(end), rec_low);
final_mod  = 100*cellfun(@(x) x(end), rec_mod);
final_high = 100*cellfun(@(x) x(end), rec_high);

final_nr   = final_nr(1:nCycles);
final_low  = final_low(1:nCycles);
final_mod  = final_mod(1:nCycles);
final_high = final_high(1:nCycles);

% ---- 2a) Line plot (new order) ----
figure; hold on;

p1 = plot(1:nCycles, final_nr,  '-sk','LineWidth',1.2,'MarkerFaceColor','k');
p2 = plot(1:nCycles, final_low, '-dg','LineWidth',1.2,'MarkerFaceColor','g');
p3 = plot(1:nCycles, final_mod, '-ob','LineWidth',1.2,'MarkerFaceColor','b');
p4 = plot(1:nCycles, final_high,'-^m','LineWidth',1.2,'MarkerFaceColor','m');

xlabel('Cycle');
ylabel('Final Recovery [%]');
title('Final H_2 Recovery per Cycle');
grid on;
legend([p1 p2 p3 p4], ...
    {'No Reaction', ...
     'Reaction - Low Rate', ...
     'Reaction - Moderate Rate', ...
     'Reaction - High Rate'}, ...
    'Location','best');

% ---- 2b) Grouped bar differences relative to No Reaction ----
diff_low  = final_low  - final_nr;
diff_mod  = final_mod  - final_nr;
diff_high = final_high - final_nr;

diffMat = [diff_low(:), diff_mod(:), diff_high(:)];

figure;
hb = bar(1:nCycles, diffMat, 'grouped');
xlabel('Cycle'); ylabel('\Delta Recovery relative to No Reaction [%]');
title('Effect of Microbial Reactions on H_2 Recovery');

grid on;

legend(hb, ...
    {'Low Rate - No Reaction', ...
     'Moderate Rate - No Reaction', ...
     'High Rate - No Reaction'}, ...
    'Location','best');



%% ============================================================
% 3) Purity vs time – all cycles (moderate case like original)
% ============================================================
figure; hold on;
for c = 1:nCycles
    plot(time_mod{c}, 100*pur_mod{c}, 'LineWidth', 1.2, ...
        'DisplayName', sprintf('Cycle %d', c));
end
xlabel('Date'); ylabel('H_2 Purity [%]');
title('H_2 Purity vs Time (Moderate Rate)');
legend show; grid on;



%% ============================================================
% 4) Average Purity vs cycle (lines + grouped bars)
%    Order: No → Low → Moderate → High
% ============================================================
avg_nr   = 100*cellfun(@mean, pur_nr);
avg_low  = 100*cellfun(@mean, pur_low);
avg_mod  = 100*cellfun(@mean, pur_mod);
avg_high = 100*cellfun(@mean, pur_high);

avg_nr   = avg_nr(1:nCycles);
avg_low  = avg_low(1:nCycles);
avg_mod  = avg_mod(1:nCycles);
avg_high = avg_high(1:nCycles);

% ---- Line plot ----
figure; hold on;

p1 = plot(1:nCycles, avg_nr,  '-sk','LineWidth',1.2,'MarkerFaceColor','k');
p2 = plot(1:nCycles, avg_low, '-dg','LineWidth',1.2,'MarkerFaceColor','g');
p3 = plot(1:nCycles, avg_mod, '-ob','LineWidth',1.2,'MarkerFaceColor','b');
p4 = plot(1:nCycles, avg_high,'-^m','LineWidth',1.2,'MarkerFaceColor','m');

xlabel('Cycle');
ylabel('Average H_2 Purity [%]');
title('Average H_2 Purity per Cycle');
grid on;
legend([p1 p2 p3 p4], ...
    {'No Reaction', ...
     'Reaction - Low Rate', ...
     'Reaction - Moderate Rate', ...
     'Reaction - High Rate'}, ...
    'Location','best');

% ---- Difference bars ----
diff_low  = avg_low  - avg_nr;
diff_mod  = avg_mod  - avg_nr;
diff_high = avg_high - avg_nr;

diffMat = [diff_low(:), diff_mod(:), diff_high(:)];

figure;
hb = bar(1:nCycles, diffMat, 'grouped');
xlabel('Cycle');
ylabel('\Delta Purity relative to No Reaction [%]');
title('Impact of Microbial Reactions on Average H_2 Purity');
grid on;

legend(hb, ...
    {'Low Rate - No Reaction', ...
     'Moderate Rate - No Reaction', ...
     'High Rate - No Reaction'}, ...
    'Location','best');



%% ============================================================
% 5) H2S vs time (moderate case)
% ============================================================
figure; hold on;
for c = 1:nCycles
    plot(time_mod{c}, H2S_mod{c}, 'LineWidth', 1.2, ...
        'DisplayName', sprintf('Cycle %d', c));
end
xlabel('Date'); ylabel('H_2S Concentration [ppm]');
title('H_2S Concentration vs Time (Moderate Rate)');
legend show; grid on;



%% ============================================================
% 6) Avg H2S vs cycle (only exists for reaction cases)
% ============================================================
figure; hold on;
plot(1:nCycles, avgH2S_mod, '-b','LineWidth',1.2);
plot(1:nCycles, avgH2S_mod, 'or','MarkerSize',6,'MarkerFaceColor','r');
xlabel('Cycle'); ylabel('Average H_2S [ppm]');
title('Average H_2S per Cycle (Moderate Rate Only)');
grid on;



%% ============================================================
% 7) BHP plot – all datasets, correct order + labels
% ============================================================
ctrl = schedule.step.control;
dCtrl = diff([0; ctrl==1]);
injStart = find(dCtrl==1);
dCtrl = diff([ctrl==3; 0]);
prodEnd = find(dCtrl==-1);

nCycles = min([numel(injStart), numel(prodEnd), ...
               numel(bhp_nr), numel(bhp_low), ...
               numel(bhp_mod), numel(bhp_high)]);

figure('Units','normalized','Position',[0.05 0.2 0.85 0.5]); hold on;

for c = 1:nCycles
    range = injStart(c):prodEnd(c);

    plot(stepDates(range), bhp_nr{c},  '-k', 'LineWidth',1.2, ...
         'DisplayName', sprintf('Cycle %d (No Reaction)', c));

    plot(stepDates(range), bhp_low{c},  '--g', 'LineWidth',1.2, ...
         'DisplayName', sprintf('Cycle %d (Low Rate)', c));

    plot(stepDates(range), bhp_mod{c},  '-.b', 'LineWidth',1.2, ...
         'DisplayName', sprintf('Cycle %d (Moderate Rate)', c));

    plot(stepDates(range), bhp_high{c}, ':m', 'LineWidth',1.4, ...
         'DisplayName', sprintf('Cycle %d (High Rate)', c));
end

xlabel('Date');
ylabel('BHP [bar]');
title('Bottomhole Pressure vs Time (All Scenarios)');
legend('Location','eastoutside');
grid on;
