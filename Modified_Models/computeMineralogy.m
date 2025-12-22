function [states, mineralData] = computeMineralogy(model, states, t)
%COMPUTEMINERALOGY  Compute mineral moles in each cell (from wt% at t=1).
%
% Inputs:
%   model  - MRST model (uses model.G, model.rock)
%   states - state cell array (uses states{t}.Mineralogy.*)
%   t      - current time step (1-based)
%
% Outputs:
%   states      - updated states (Mineralogy.* at step t converted to moles for t=1)
%   mineralData - struct with per–cell rock volume, density, total mass, and
%                 initial mineral moles (Calcite_i, Dolomite_i, Quartz_i, ...)

    G    = model.G;
    rock = model.rock;
    nc   = G.cells.num;

    % Preallocate
    V_Rock        = zeros(nc,1);   % rock volume per cell [m^3]
    rho_ave       = zeros(nc,1);   % average rock density [g/cm^3]
    m_tot         = zeros(nc,1);   % total rock mass per cell [g]

    Quartz_i      = zeros(nc,1);
    Calcite_i     = zeros(nc,1);
    Dolomite_i    = zeros(nc,1);
    Anhydrite_i   = zeros(nc,1);
    Goethite_i    = zeros(nc,1);
    Brucite_i     = zeros(nc,1);
    Portlandite_i = zeros(nc,1);
    Pyrite_i      = zeros(nc,1);
    Gypsum_i      = zeros(nc,1);

    for i = 1:nc
        % Rock volume in cell i
        V_Rock(i) = (1 - rock.poro(i)) * G.cells.volumes(i);   % [m^3]

        if t == 1
            % Weight fractions from Mineralogy (given in wt%)
            Quartz_Share      = states{t,1}.Mineralogy.Quartz(i)      * 0.01;
            Calcite_Share     = states{t,1}.Mineralogy.Calcite(i)     * 0.01;
            Dolomite_Share    = states{t,1}.Mineralogy.Dolomite(i)    * 0.01;
            Anhydrite_Share   = states{t,1}.Mineralogy.Anhydrite(i)   * 0.01;
            Goethite_Share    = states{t,1}.Mineralogy.Goethite(i)    * 0.01;
            Brucite_Share     = states{t,1}.Mineralogy.Brucite(i)     * 0.01;
            Portlandite_Share = states{t,1}.Mineralogy.Portlandite(i) * 0.01;
            Pyrite_Share      = states{t,1}.Mineralogy.Pyrite(i)      * 0.01;
            Gypsum_Share      = states{t,1}.Mineralogy.Gypsum(i)      * 0.01;

            % Total weight fraction (should be ~1)
            A = Quartz_Share + Calcite_Share + Dolomite_Share + ...
                Anhydrite_Share + Goethite_Share + Brucite_Share + ...
                Portlandite_Share + Pyrite_Share + Gypsum_Share;

            % Weighted by specific densities [g/cm^3]
            % (kept identical to your original B expression)
            B = Quartz_Share   / 2.65  + ...
                Calcite_Share  / 2.71  + ...
                Dolomite_Share / 2.84  + ...
                Anhydrite_Share/ 2.97  + ...
                Goethite_Share / 4.26  + ...
                Brucite_Share  / 2.40  + ...
                Pyrite_Share   / 4.9   + ...
                Gypsum_Share   / 2.39;

            rho_ave(i) = A / B;              % [g/cm^3]
            m_tot(i)   = rho_ave(i) * V_Rock(i) * 1000; % [g] (1 m^3 = 1e6 cm^3)

            % Convert mass [g] to moles per cell, then to mmol (x1000)
            Quartz_i(i)      = Quartz_Share      * m_tot(i) / 60.083   * 1000;
            Calcite_i(i)     = Calcite_Share     * m_tot(i) / 100.0869 * 1000;
            Dolomite_i(i)    = Dolomite_Share    * m_tot(i) / 184.40   * 1000;
            Anhydrite_i(i)   = Anhydrite_Share   * m_tot(i) / 136.14   * 1000;
            Goethite_i(i)    = Goethite_Share    * m_tot(i) / 88.85    * 1000;
            Brucite_i(i)     = Brucite_Share     * m_tot(i) / 58.32    * 1000;
            Portlandite_i(i) = Portlandite_Share * m_tot(i) / 74.09    * 1000;
            Pyrite_i(i)      = Pyrite_Share      * m_tot(i) / 119.98   * 1000;
            Gypsum_i(i)      = Gypsum_Share      * m_tot(i) / 172.2    * 1000;

            % Overwrite Mineralogy at t=1 with moles (what your original code does)
            states{t,1}.Mineralogy.Calcite(i)     = Calcite_i(i);
            states{t,1}.Mineralogy.Dolomite(i)    = Dolomite_i(i);
            states{t,1}.Mineralogy.Quartz(i)      = Quartz_i(i);
            states{t,1}.Mineralogy.Anhydrite(i)   = Anhydrite_i(i);
            states{t,1}.Mineralogy.Goethite(i)    = Goethite_i(i);
            states{t,1}.Mineralogy.Brucite(i)     = Brucite_i(i);
            states{t,1}.Mineralogy.Portlandite(i) = Portlandite_i(i);
            states{t,1}.Mineralogy.Pyrite(i)      = Pyrite_i(i);
            states{t,1}.Mineralogy.Gypsum(i)      = Gypsum_i(i);

        else
            % For t > 1 use previous step moles as the "initial" values
            Calcite_i(i)     = states{t-1,1}.Mineralogy.Calcite(i);
            Dolomite_i(i)    = states{t-1,1}.Mineralogy.Dolomite(i);
            Quartz_i(i)      = states{t-1,1}.Mineralogy.Quartz(i);
            Anhydrite_i(i)   = states{t-1,1}.Mineralogy.Anhydrite(i);
            Goethite_i(i)    = states{t-1,1}.Mineralogy.Goethite(i);
            Brucite_i(i)     = states{t-1,1}.Mineralogy.Brucite(i);
            Portlandite_i(i) = states{t-1,1}.Mineralogy.Portlandite(i);
            Pyrite_i(i)      = states{t-1,1}.Mineralogy.Pyrite(i);
            Gypsum_i(i)      = states{t-1,1}.Mineralogy.Gypsum(i);
        end
    end

    % Pack outputs
    mineralData = struct( ...
        'V_Rock',        V_Rock, ...
        'rho_ave',       rho_ave, ...
        'm_tot',         m_tot, ...
        'Quartz_i',      Quartz_i, ...
        'Calcite_i',     Calcite_i, ...
        'Dolomite_i',    Dolomite_i, ...
        'Anhydrite_i',   Anhydrite_i, ...
        'Goethite_i',    Goethite_i, ...
        'Brucite_i',     Brucite_i, ...
        'Portlandite_i', Portlandite_i, ...
        'Pyrite_i',      Pyrite_i, ...
        'Gypsum_i',      Gypsum_i ...
    );
end
