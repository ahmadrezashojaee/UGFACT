function [WaterTable, GasTable] = coreyBrooksTables(Swi, Sgc, krw_end, krg_end, nw, ng, nPoints)
% COREYBROOKSTABLES: Generate Corey–Brooks rel perm tables for gas-water system.
%
% Inputs:
%   Swi      - irreducible water saturation
%   Sgc      - critical/residual gas saturation
%   krw_end  - endpoint rel perm for water
%   krg_end  - endpoint rel perm for gas
%   nw, ng   - Corey exponents (water, gas)
%   nPoints  - number of points in the table
%
% Outputs:
%   WaterTable - [Sw, krw]
%   GasTable   - [Sg, krg]

    % Water saturation
    Sw = linspace(Swi, 1 - Sgc, nPoints)';
    Sg = 1 - Sw;

    % Effective water saturation
    Swe = (Sw - Swi) ./ max(1e-12, (1 - Swi - Sgc));
    Swe = min(max(Swe, 0), 1);

    % Corey-Brooks rel perms
    krw = krw_end .* (Swe .^ nw);
    krg = krg_end .* ((1 - Swe) .^ ng);

    % Build tables
    WaterTable_main = [Sw, krw];
    GasTable_main   = [Sg, krg];

    % --- Add extra rows at Sw/Sg = 0 and 1 ---
    WaterTable = [
        0,            0;          % Sw = 0, krw = 0
        WaterTable_main;
        1,            krw_end     % Sw = 1, krw = endpoint
    ];

    GasTable = [
        1,            krg_end          % Sg = 1, krg = endpoint
        GasTable_main;
        0,            0                % Sg = 1, krg = 0
    ];
end
