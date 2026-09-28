function addBigColorbar(ax)
    % Make room for a wider colorbar
    ax.Units = 'normalized';
    axPos = ax.Position;

    % Slightly shrink axes width to create space
    ax.Position = [axPos(1), axPos(2), axPos(3)*0.86, axPos(4)];

    % Create colorbar
    cb = colorbar(ax, 'eastoutside');
    cb.FontSize = 10;
    cb.Units = 'normalized';

    % Place colorbar manually
    axPos = ax.Position;
    cb.Position = [axPos(1)+axPos(3)+0.006, ...   % x
                   axPos(2)+0.03*axPos(4), ...    % y
                   0.012, ...                      % width
                   0.94*axPos(4)];                % height
end