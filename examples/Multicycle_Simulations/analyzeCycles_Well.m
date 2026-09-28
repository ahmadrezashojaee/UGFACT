function [purityCycles, recoveryCycles, timeCycles, H2SCycles, avgH2S, bhpCycles, injH2Cycles] = analyzeCycles_Well(wellSol, schedule, stepDates)
% analyzeCycles_Well
% ---------------------------------------------------------
% Compute H2 purity, recovery, H2S concentration, BHP, and injected H2 per cycle
%
% INPUTS:
%   wellSol   - cell array of well solutions {t}
%   schedule  - MRST schedule struct
%   stepDates - datetime array [nSteps,1]
%
% OUTPUTS:
%   purityCycles   - cell array, each entry vector of H2 purity over prod steps
%   recoveryCycles - cell array, each entry vector of recovery over prod steps
%   timeCycles     - cell array of datetimes for plotting
%   H2SCycles      - cell array, each entry vector of H2S concentration (ppm)
%   avgH2S         - vector, average H2S (ppm) per cycle
%   bhpCycles      - cell array, each entry vector of BHP (bar) over whole cycle
%   injH2Cycles    - vector, total injected H2 (mol) per cycle
%
% Note:
%   Injection phase = ctrl == 1
%   Production phase = ctrl == 3

% --- Molecular weights (kg/mol) ---
MW.H2  = 2.016e-3; 
MW.CH4 = 16.04e-3; 
MW.CO2 = 44.01e-3; 
MW.H2S = 34.08e-3;

ctrl = schedule.step.control;

% --- Detect injection starts and production ends ---
dCtrl = diff([0; ctrl==1]);
injStart = find(dCtrl==1);        % start index of injection
dCtrl = diff([ctrl==3; 0]);
prodEnd = find(dCtrl==-1);        % end index of production

nCycles = min(numel(injStart), numel(prodEnd));

purityCycles   = cell(nCycles,1);
recoveryCycles = cell(nCycles,1);
timeCycles     = cell(nCycles,1);
H2SCycles      = cell(nCycles,1);
avgH2S         = zeros(nCycles,1);
bhpCycles      = cell(nCycles,1);
injH2Cycles    = zeros(nCycles,1); % NEW: store injected H2 per cycle

for c = 1:nCycles
    range = injStart(c):prodEnd(c);

    % --- Injection & Production indices ---
    idxInj  = range(ctrl(range)==1);
    idxProd = range(ctrl(range)==3);

    % ---- Total injected H2 (mol) ----
    injH2 = 0;
    for t = idxInj(:)'
        rateH2 = wellSol{t}.H2 ./ MW.H2; % mol/s
        if t > 1
            dtSec = seconds(stepDates(t) - stepDates(t-1));
        else
            dtSec = seconds(stepDates(t) - stepDates(1));
        end
        injH2 = injH2 + rateH2*dtSec;
    end
    injH2Cycles(c) = injH2; % save per cycle

    % ---- Production phase ----
    cumProdH2 = 0;
    cumProdAll = 0;
    purity = zeros(numel(idxProd),1);
    recovery = zeros(numel(idxProd),1);
    H2Sppm = zeros(numel(idxProd),1);

    for k = 1:numel(idxProd)
        t = idxProd(k);

        % Production rates (mol/s, negative sign convention)
        qH2  = -wellSol{t}.H2  ./ MW.H2;
        qCH4 = -wellSol{t}.CH4 ./ MW.CH4;
        qCO2 = -wellSol{t}.CO2 ./ MW.CO2;
        qH2S = -wellSol{t}.H2S ./ MW.H2S;

        qTot = qH2 + qCH4 + qCO2 + qH2S;

        if t > 1
            dtSec = seconds(stepDates(t) - stepDates(t-1));
        else
            dtSec = seconds(stepDates(t) - stepDates(1));
        end

        % Integrate
        cumProdH2  = cumProdH2 + qH2*dtSec;
        cumProdAll = cumProdAll + qTot*dtSec;

        % Purity, recovery, H2S concentration
        purity(k)   = qH2 / max(qTot,eps);
        recovery(k) = cumProdH2 / max(injH2,eps);
        H2Sppm(k)   = (qH2S / max(qTot,eps)) * 1e6; % ppm
    end

    % --- Save per-cycle outputs ---
    purityCycles{c}   = purity;
    recoveryCycles{c} = recovery;
    timeCycles{c}     = stepDates(idxProd);
    H2SCycles{c}      = H2Sppm;
    avgH2S(c)         = mean(H2Sppm);

    % --- Extract BHP over whole cycle (inj + prod) ---
    bhpBar = zeros(numel(range),1);
    for k = 1:numel(range)
        t = range(k);
        bhpBar(k) = wellSol{t}.bhp / 1e5; % Pa → bar
    end
    bhpCycles{c} = bhpBar;
end
end
