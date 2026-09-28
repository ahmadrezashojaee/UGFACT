% Define vectors
Sw_vec = [30, 45, 60];
Mu_vec = unique([0, 0.3:0.5:4.1 4.1]);  % Ensure no duplicate values
time = (1:200)';

% Generate turbo colormap (blue→red)
colors = turbo(length(Mu_vec));  % [N×3] RGB array

% Preallocate array to hold the ΔP axes handles
ax2 = gobjects(length(Sw_vec),1);

% Initialize Excel file
excel_filename = 'pressure_analysis.xlsx';
if isfile(excel_filename)
    delete(excel_filename);
end

% Create figure
figure;
for i = 1:length(Sw_vec)
    Sw_val = Sw_vec(i);

    % Containers for pressure data
    pressure_map = containers.Map('KeyType','char','ValueType','any');
    base_p = [];
    legendEntries_dp = {};

    % Tables for Excel export
    pressure_table = table(time, 'VariableNames', {'Time'});
    deltaP_table    = table(time, 'VariableNames', {'Time'});

    % Load all pressures
    for j = 1:length(Mu_vec)
        Mu_val = Mu_vec(j);
        % format Mu string
        if Mu_val == 0
            Mu_str = '0';
        else
            Mu_str = num2str(Mu_val,'%.1f');
        end

        % construct filename
        res_filename = sprintf('Sw_%d_Mu_%s.mat', Sw_val, Mu_str);
        if ~isfile(res_filename)
            warning('File %s not found.', res_filename);
            continue;
        end

        % load data
        S = load(res_filename);
        varName = fieldnames(S);
        resData = S.(varName{1});

        % extract pressure
        p = zeros(length(time),1);
        for t = 1:length(time)
            p(t) = resData{t,1}.pressure(1)/1e5;
        end

        % store
        pressure_map(Mu_str) = p;
        pressure_table.(sprintf('mu_max_%.1f',Mu_val)) = p;
        if Mu_val==0
            base_p = p;
        end
    end

    if isempty(base_p)
        warning('Missing base mu_max = 0 for Sw = %d. Skipping.', Sw_val);
        continue;
    end

    % --- Top row: Pressure vs time ---
    subplot(2,3,i); hold on;
    for j = 1:length(Mu_vec)
        Mu_val = Mu_vec(j);
        if Mu_val==0
            Mu_str='0'; ls='--'; lw=3;
        else
            Mu_str=num2str(Mu_val,'%.1f'); ls='-'; lw=2;
        end
        if isKey(pressure_map,Mu_str)
            plot(time, pressure_map(Mu_str), ls, 'LineWidth',lw, 'Color', colors(j,:));
        end
    end
    xlabel('Time (day)');
    ylabel('Pressure (bar)');
    title(sprintf('Sw = %d: Pressure vs Time',Sw_val));
    grid on;

    % --- Bottom row: ΔP vs time ---
    ax2(i) = subplot(2,3,i+3); hold on;
    % baseline ΔP = 0
    plot(time, zeros(size(time)), '--','LineWidth',3,'Color',colors(1,:));
    deltaP_table.('Delta_mu_max_0') = zeros(size(time));
    legendEntries_dp{1} = '\mu_{max}=0';

    for j = 2:length(Mu_vec)
        Mu_val = Mu_vec(j);
        Mu_str = num2str(Mu_val,'%.1f');
        if isKey(pressure_map,Mu_str)
            delta_p = base_p - pressure_map(Mu_str);
            plot(time, delta_p, '-','LineWidth',2,'Color',colors(j,:));
            legendEntries_dp{end+1} = sprintf('\\mu_{max}=%.1f',Mu_val);
            deltaP_table.(sprintf('Delta_mu_max_%.1f',Mu_val)) = delta_p;
        end
    end
    xlabel('Time (day)');
    ylabel('\DeltaP (bar)');
    title(sprintf('Sw = %d: \\DeltaP vs Time',Sw_val));
    grid on;

    % show legend only in last ΔP subplot
    if i==3
        legend(legendEntries_dp,'Location','best','Orientation','Horizontal');
    end

    % write to Excel
    writetable(pressure_table, excel_filename, 'Sheet', sprintf('Sw_%d_Pressure',Sw_val));
    writetable(deltaP_table,   excel_filename, 'Sheet', sprintf('Sw_%d_DeltaP',Sw_val));
end

% link all ΔP subplots to share the same y-limits
linkaxes(ax2,'y');
