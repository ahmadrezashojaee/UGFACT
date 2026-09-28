%% =================== Postprocessing: H2 Recovery & Purity ===================
% Assumes you have wellSol from MRST simulation

% ---------------- User inputs ----------------
nCycles = 20; 
nInj    = 31; 
nShut   = 13; 
nProd   = 31; 
stepsPerCycle = nInj + nShut + nProd;

dt = 5*day;          % timestep length [s]
dt_days = dt/day;    % timestep length [days]

% Molecular weights [kg/mol]
MW_H2  = 2.016e-3;
MW_CH4 = 16.04e-3;
MW_CO2 = 44.01e-3;
MW_H2S = 34.08e-3;

% ---------------- Preallocate ----------------
H2_inj_total  = zeros(nCycles,1);
H2_prod_total = zeros(nCycles,1);
Recovery      = zeros(nCycles,1);
Purity_cycle  = zeros(nCycles,1);

Purity_time    = nan(nProd, nCycles);   % time-resolved purity
H2rate_time    = nan(nProd, nCycles);   % H2 rate [mol/day]
CH4rate_time   = nan(nProd, nCycles);   % CH4 rate [mol/day]
time_axis      = (1:nProd)*dt_days;     % days within a production stage

% ---------------- Loop over cycles ----------------
for c = 1:nCycles
    injSteps  = (c-1)*stepsPerCycle + (1:nInj);
    prodSteps = (c-1)*stepsPerCycle + (nInj+nShut+1 : stepsPerCycle);

    % --- Injected H2 (mol) ---
    injH2 = 0;
    for t = injSteps
        % kg/s -> mol/s -> mol per step
        injH2 = injH2 + (wellSol{t}.H2 / MW_H2) * dt;
    end
    H2_inj_total(c) = injH2;

    % --- Produced gases (mol) ---
    prodH2 = 0; 
    prodTotal = 0;
    for k = 1:nProd
        t = prodSteps(k);

        % kg/s -> mol/s -> mol per step
        h2   = (-wellSol{t}.H2  / MW_H2) * dt;
        ch4  = (-wellSol{t}.CH4 / MW_CH4) * dt;
        co2  = (-wellSol{t}.CO2 / MW_CO2) * dt;
        h2s  = (-wellSol{t}.H2S / MW_H2S) * dt;

        stepTot = max(h2+ch4+co2+h2s, 0);

        prodH2    = prodH2    + max(h2,0);
        prodTotal = prodTotal + stepTot;

        % time-resolved purity (dimensionless)
        Purity_time(k,c) = h2 / max(stepTot,1e-12);

        % ---- NEW: H2 and CH4 production rates in mol/day ----
        H2rate_time(k,c)  = max(h2,0)  / dt_days;   % [mol/day]
        CH4rate_time(k,c) = max(ch4,0) / dt_days;   % [mol/day]
    end
    H2_prod_total(c) = prodH2;

    % --- Metrics per cycle ---
    Recovery(c)     = prodH2 / injH2;
    Purity_cycle(c) = prodH2 / prodTotal;
end

%% =================== Plotting ===================

% ---- Cycle-averaged metrics ----
figure;
subplot(2,1,1);
plot(1:nCycles, Recovery*100, 'o-b','LineWidth',2);
xlabel('Cycle'); ylabel('H2 Recovery (%)');
grid on; title('Cycle-wise H2 Recovery');

subplot(2,1,2);
plot(1:nCycles, Purity_cycle*100, 's-r','LineWidth',2);
xlabel('Cycle'); ylabel('H2 Purity (%)');
grid on; title('Cycle-wise H2 Purity');

% ---- Time-resolved purity per cycle ----
figure; hold on;
for c = 1:nCycles
    plot(time_axis, Purity_time(:,c)*100, 'LineWidth',1.5, ...
         'DisplayName', sprintf('Cycle %d',c));
end
xlabel('Production time (days)'); ylabel('H2 Purity (%)');
title('H2 Purity vs Time for Each Cycle');
legend show; grid on;

% ---- Continuous timeline (purity, all cycles) ----
figure; hold on;
for c = 1:nCycles
    t0 = (c-1)*stepsPerCycle*dt_days;  % start time of this cycle [days]
    plot(t0 + time_axis, Purity_time(:,c)*100, 'LineWidth',1.5);
end
xlabel('Simulation time (days)'); ylabel('H2 Purity (%)');
title('H2 Purity over Simulation'); grid on;

%% ===== H2 molar production rate vs time (mol/day) =====

% ---- H2 rate vs production time for each cycle ----
figure; hold on;
for c = 1:nCycles
    plot(time_axis, H2rate_time(:,c), 'LineWidth',1.5, ...
         'DisplayName', sprintf('Cycle %d', c));
end
xlabel('Production time (days)');
ylabel('H2 production rate (mol/day)');
title('H2 production rate vs time for each cycle');
legend show; grid on;

% ---- H2 rate over full simulation (all cycles concatenated) ----
figure; hold on;
for c = 1:nCycles
    t0 = (c-1)*stepsPerCycle*dt_days;  % start time of this cycle [days]
    plot(t0 + time_axis, H2rate_time(:,c), 'LineWidth',1.5);
end
xlabel('Simulation time (days)');
ylabel('H2 production rate (mol/day)');
title('H2 production rate over full simulation');
grid on;

%% ===== Export H2 and CH4 rates (mol/day) to Excel (with REAL DATES) =====

RealStart = datetime(2025,1,1);   % choose real start date

totalProdSteps = nCycles * nProd;
Time_full      = nan(totalProdSteps,1);   % linear days
H2_full        = nan(totalProdSteps,1);
CH4_full       = nan(totalProdSteps,1);
RealDate_full  = NaT(totalProdSteps,1);   % real calendar time

idx = 0;
for c = 1:nCycles
    t0 = (c-1)*stepsPerCycle*dt_days;  % start time of this cycle [days]
    for k = 1:nProd
        idx = idx + 1;

        % Linear simulation time (days)
        Time_full(idx) = t0 + time_axis(k);

        % Real calendar date
        RealDate_full(idx) = RealStart + days(Time_full(idx));

        % Save rates
        H2_full(idx)  = H2rate_time(k,c);
        CH4_full(idx) = CH4rate_time(k,c);
    end
end

% Create table
T = table(RealDate_full, Time_full, H2_full, CH4_full, ...
          'VariableNames', {'RealDate','Time_days','H2_mol_per_day','CH4_mol_per_day'});

% Write to Excel
writetable(T, 'H2_CH4_rates_with_real_dates.xlsx');
