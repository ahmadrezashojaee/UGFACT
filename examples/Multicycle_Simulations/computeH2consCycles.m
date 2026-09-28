function H2consCycles = computeH2consCycles(statesR, stepDates, nPerCycle)

    NC      = size(statesR{1}.y,1);
    nSteps  = numel(statesR);
    nCycles = nSteps / nPerCycle;

    H2consCycles = zeros(NC, nCycles);   % mol per cell per cycle

    for c = 1:nCycles
        idxRange = (c-1)*nPerCycle + (1:nPerCycle);  % steps in this cycle
        H2cycle  = zeros(NC,1);                     % mol in this cycle

        for k = 2:numel(idxRange)
            t     = idxRange(k);
            tprev = idxRange(k-1);

            dt = days(stepDates(t) - stepDates(tprev));   % days
            if dt <= 0
                continue;
            end

            st = statesR{t};
            if ~isfield(st,'Solution')
                continue;
            end

            % Rates [mmol/day/kg_w]
            r_MET = st.Solution.MET_Rate;
            r_SRB = st.Solution.SRB_Rate;
            r_ACE = st.Solution.ACE_Rate;

            r_tot = r_MET + r_SRB + r_ACE;   % total mmol/day/kg_w

            % Water mass per cell [kg]
            m_water = st.FlowProps.ComponentTotalMass{1,1};

            % H2 consumed this step: r_tot * m_water * dt
            % [mmol/day/kg * kg * day = mmol] -> mol
            H2cycle = H2cycle + r_tot .* m_water * dt * 1e-3;
        end

        H2consCycles(:,c) = H2cycle;
    end
end
