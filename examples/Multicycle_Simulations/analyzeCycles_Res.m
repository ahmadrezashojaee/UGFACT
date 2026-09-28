function [H2consCycles, timeCycles, H2consStep, H2consTotal, others] = analyzeCycles_Res(states, schedule, stepDates)
% analyzeCycles_Res
% ---------------------------------------------------------
% Computes hydrogen consumption per cycle (reset each cycle)
% and global cumulative totals across the whole simulation.
% Also computes other useful diagnostics (e.g. average pressure).
%
% INPUTS:
%   states   - cell array with reservoir states {t,1}
%   schedule - MRST schedule struct
%   stepDates- vector of datetimes for simulation timesteps
%
% OUTPUTS:
%   H2consCycles - cell array (nCycles x 1), cumulative H2 consumption per cycle
%   timeCycles   - cell array (nCycles x 1), time in days per cycle
%   H2consStep   - cell array (nSteps x 1), per-grid struct
%   H2consTotal  - struct with fields (global cumulative across simulation)
%   others       - struct with additional diagnostics:
%                  .Average_Pres : struct with
%                     .Values -> vector of average pressure [bar] per step
%                     .Time   -> datetime vector of step dates

    nSteps = numel(states);

    % --- Identify cycle boundaries from schedule ---
    ctrl = schedule.step.control;
    cycleEndIdx = find(ctrl(1:end-1) == 3 & ctrl(2:end) == 1); % 3->1 transition
    cycleEndIdx = [cycleEndIdx; nSteps]; % last cycle goes to end
    cycleStartIdx = [1; cycleEndIdx(1:end-1)+1];
    nCycles = numel(cycleStartIdx);

    % --- Outputs ---
    H2consCycles = cell(nCycles,1);
    timeCycles   = cell(nCycles,1);
    H2consStep   = cell(nSteps,1);

    % global reservoir totals (cumulative across simulation)
    MET_tot  = zeros(nSteps,1);
    SRB_tot  = zeros(nSteps,1);
    ACE_tot  = zeros(nSteps,1);
    TOT_tot  = zeros(nSteps,1);

    % --- NEW: Average pressure vector ---
    avgPresVals = zeros(nSteps,1);

    % --- Loop over all steps ---
    for t = 1:nSteps
        if t == 1
            dt = 0;
        else
            dt = days(stepDates(t) - stepDates(t-1));
        end

        % rates [mmol/day/kgw]
        if isfield(states{1,1}, 'Solution')
            r_MET = states{t}.Solution.MET_Rate;
            r_SRB = states{t}.Solution.SRB_Rate;
            r_ACE = states{t}.Solution.ACE_Rate;
        else
            r_MET = 0;
            r_SRB = 0;
            r_ACE = 0;
        end

        % water mass [kg]
        m_water = states{t}.FlowProps.ComponentTotalMass{1,1};

        % per grid block consumption [mol]
        H2_MET = r_MET .* m_water * dt / 1000;
        H2_SRB = r_SRB .* m_water * dt / 1000;
        H2_ACE = r_ACE .* m_water * dt / 1000;
        H2_TOT = H2_MET + H2_SRB + H2_ACE;

        % store per-grid results
        H2consStep{t} = struct( ...
            'MET',  H2_MET(:), ...
            'SRB',  H2_SRB(:), ...
            'ACE',  H2_ACE(:), ...
            'Total',H2_TOT(:));

        % global cumulative totals
        MET_step = sum(H2_MET);
        SRB_step = sum(H2_SRB);
        ACE_step = sum(H2_ACE);
        TOT_step = MET_step + SRB_step + ACE_step;

        if t == 1
            MET_tot(t) = MET_step;
            SRB_tot(t) = SRB_step;
            ACE_tot(t) = ACE_step;
            TOT_tot(t) = TOT_step;
        else
            MET_tot(t) = MET_tot(t-1) + MET_step;
            SRB_tot(t) = SRB_tot(t-1) + SRB_step;
            ACE_tot(t) = ACE_tot(t-1) + ACE_step;
            TOT_tot(t) = TOT_tot(t-1) + TOT_step;
        end

        % --- NEW: Average reservoir pressure (Pa → bar) ---
        avgPresVals(t) = mean(states{t}.pressure(:)) / 1e5;
    end

    % --- Cycle-wise cumulative (reset each cycle) ---
    for c = 1:nCycles
        idxRange = cycleStartIdx(c):cycleEndIdx(c);
        nLocal   = numel(idxRange);

        MET_cum = zeros(nLocal,1);
        SRB_cum = zeros(nLocal,1);
        ACE_cum = zeros(nLocal,1);
        TOT_cum = zeros(nLocal,1);
        timeLoc = zeros(nLocal,1);

        for k = 1:nLocal
            t = idxRange(k);

            MET_step = sum(H2consStep{t}.MET);
            SRB_step = sum(H2consStep{t}.SRB);
            ACE_step = sum(H2consStep{t}.ACE);
            TOT_step = sum(H2consStep{t}.Total);

            if k == 1
                MET_cum(k) = MET_step;
                SRB_cum(k) = SRB_step;
                ACE_cum(k) = ACE_step;
                TOT_cum(k) = TOT_step;
                timeLoc(k) = 0;
            else
                MET_cum(k) = MET_cum(k-1) + MET_step;
                SRB_cum(k) = SRB_cum(k-1) + SRB_step;
                ACE_cum(k) = ACE_cum(k-1) + ACE_step;
                TOT_cum(k) = TOT_cum(k-1) + TOT_step;
                timeLoc(k) = days(stepDates(t) - stepDates(idxRange(1)));
            end
        end

        H2consCycles{c} = struct( ...
            'MET',  MET_cum, ...
            'SRB',  SRB_cum, ...
            'ACE',  ACE_cum, ...
            'Total',TOT_cum);

        timeCycles{c} = timeLoc;
    end

    % --- Global cumulative totals ---
    H2consTotal = struct( ...
        'MET',  MET_tot, ...
        'SRB',  SRB_tot, ...
        'ACE',  ACE_tot, ...
        'Total',TOT_tot);

    % --- NEW: Others output ---
    others = struct();
    others.Average_Pres = struct( ...
        'Values', avgPresVals, ...
        'Time',   stepDates );
end
