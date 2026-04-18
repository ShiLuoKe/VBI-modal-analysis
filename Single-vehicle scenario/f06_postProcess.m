%% f06_postProcess.m
% =========================================================================
% Description : Post-processing and visualization of the VBI system results.
%               Generates four figures:
%               Fig.1 - Vehicle & bridge mid-span displacement vs time
%               Fig.2 - Vehicle & bridge mid-span velocity vs time
%               Fig.3 - Vehicle & bridge mid-span acceleration vs time
%               Fig.4 - Coupled bridge frequency and damping ratio vs time
% Author      : Shi LuoKe (史罗克)
% Date        : 2026-04-08
% Acknowledgment: To improve code readability, this code was organized
%                 and formatted with the assistance of Claude (Anthropic).
%                 The correctness of the code has been verified by the author.
% -------------------------------------------------------------------------
% Input  : eigenResult - Coupled eigenvalue results (from subfunction f04)
%          dynResult   - Dynamic response results   (from subfunction f05)
%          theorFreq   - Theoretical uncoupled frequencies (from subfunction f02)
%          sysParams   - System parameter structure (from subfunction f01 & f02)
% =========================================================================

function f06_postProcess(eigenResult, dynResult, theorFreq, sysParams)

    %% Unpack
    dt        = sysParams.solver.dt;
    n         = sysParams.bridge.n;
    gamma_val = sysParams.solver.gamma;   % auto-detect gamma value
    Num_steps = length(dynResult.vehicle.d);
    t_dyn     = (0:Num_steps-1) * dt;

    % Theoretical reference values
    fb_ref  = theorFreq.bridge.f1;      % Uncoupled bridge frequency   [Hz]
    xib_ref = sysParams.bridge.xib;     % Uncoupled bridge damping ratio  [-]

    % ---- Common style settings ----
    CQUBlue   = '#02529f';
    fontName  = 'Times New Roman';
    fontSize  = 11;
    lineW     = 1;
    figW      = 14;   % figure width  [cm]
    figH      = 6;    % figure height [cm]

    % Legend string — auto from gamma
    legendStr = sprintf('$\\gamma = %g$', gamma_val);

    % Marker indices — equally spaced, 8 markers
    numMarkers = 8;
    idxDyn = round(linspace(1, Num_steps, numMarkers));
    idxEig = round(linspace(1, length(eigenResult.t), numMarkers));

    %% ================================================================
    %  Figure 1: Displacement — vehicle (left) & bridge mid-span (right)
    %% ================================================================
    figure('Name', 'Displacement', 'NumberTitle', 'off');
    tiledlayout(1, 2, 'Padding', 'tight', 'TileSpacing', 'tight');

    % --- (a) Vehicle displacement (hollow marker) ---
    nexttile(1);
    hold on;
    l1 = plot(t_dyn, dynResult.vehicle.d, '-s', ...
        'MarkerIndices', idxDyn, 'MarkerFaceColor', 'white', ...
        'Color', CQUBlue, 'LineWidth', lineW);
    lg = legend(l1, legendStr, ...
        'Interpreter', 'latex', 'Location', 'north', ...
        'NumColumns', 1, 'Box', 'off');
    lg.ItemTokenSize = [10, 6];
    xlabel('Time (s)', 'FontName', fontName, 'FontSize', fontSize);
    ylabel('Vehicle vertical displacement (m)', 'FontName', fontName, 'FontSize', fontSize);
    box on; grid off;
    set(gca, 'FontName', fontName, 'FontSize', fontSize);

    % --- (b) Bridge midspan displacement (filled marker) ---
    nexttile(2);
    hold on;
    l1 = plot(t_dyn, dynResult.bridge.d, '-s', ...
        'MarkerIndices', idxDyn, 'MarkerFaceColor', CQUBlue, ...
        'Color', CQUBlue, 'LineWidth', lineW);
    lg = legend(l1, legendStr, ...
        'Interpreter', 'latex', 'Location', 'north', ...
        'NumColumns', 1, 'Box', 'off');
    lg.ItemTokenSize = [10, 6];
    xlabel('Time (s)', 'FontName', fontName, 'FontSize', fontSize);
    ylabel('Bridge midspan displacement (m)', 'FontName', fontName, 'FontSize', fontSize);
    box on; grid off;
    set(gca, 'FontName', fontName, 'FontSize', fontSize);

    % ---- Figure size ----
    set(gcf, 'Units', 'centimeters', 'Position', [2 2 figW figH]);
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', fontSize);

    %% ================================================================
    %  Figure 2: Velocity — vehicle (left) & bridge mid-span (right)
    %% ================================================================
    figure('Name', 'Velocity', 'NumberTitle', 'off');
    tiledlayout(1, 2, 'Padding', 'tight', 'TileSpacing', 'tight');

    % --- (a) Vehicle velocity (hollow marker) ---
    nexttile(1);
    hold on;
    l1 = plot(t_dyn, dynResult.vehicle.v, '-s', ...
        'MarkerIndices', idxDyn, 'MarkerFaceColor', 'white', ...
        'Color', CQUBlue, 'LineWidth', lineW);
    lg = legend(l1, legendStr, ...
        'Interpreter', 'latex', 'Location', 'north', ...
        'NumColumns', 1, 'Box', 'off');
    lg.ItemTokenSize = [10, 6];
    xlabel('Time (s)', 'FontName', fontName, 'FontSize', fontSize);
    ylabel('Vehicle vertical velocity (m/s)', 'FontName', fontName, 'FontSize', fontSize);
    box on; grid off;
    set(gca, 'FontName', fontName, 'FontSize', fontSize);

    % --- (b) Bridge midspan velocity (filled marker) ---
    nexttile(2);
    hold on;
    l1 = plot(t_dyn, dynResult.bridge.v, '-s', ...
        'MarkerIndices', idxDyn, 'MarkerFaceColor', CQUBlue, ...
        'Color', CQUBlue, 'LineWidth', lineW);
    lg = legend(l1, legendStr, ...
        'Interpreter', 'latex', 'Location', 'north', ...
        'NumColumns', 1, 'Box', 'off');
    lg.ItemTokenSize = [10, 6];
    xlabel('Time (s)', 'FontName', fontName, 'FontSize', fontSize);
    ylabel('Bridge midspan velocity (m/s)', 'FontName', fontName, 'FontSize', fontSize);
    box on; grid off;
    set(gca, 'FontName', fontName, 'FontSize', fontSize);

    % ---- Figure size ----
    set(gcf, 'Units', 'centimeters', 'Position', [2 2 figW figH]);
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', fontSize);

    %% ================================================================
    %  Figure 3: Acceleration — vehicle (left) & bridge mid-span (right)
    %% ================================================================
    figure('Name', 'Acceleration', 'NumberTitle', 'off');
    tiledlayout(1, 2, 'Padding', 'tight', 'TileSpacing', 'tight');

    % --- (a) Vehicle acceleration (hollow marker) ---
    nexttile(1);
    hold on;
    l1 = plot(t_dyn, dynResult.vehicle.a, '-s', ...
        'MarkerIndices', idxDyn, 'MarkerFaceColor', 'white', ...
        'Color', CQUBlue, 'LineWidth', lineW);
    lg = legend(l1, legendStr, ...
        'Interpreter', 'latex', 'Location', 'north', ...
        'NumColumns', 1, 'Box', 'off');
    lg.ItemTokenSize = [10, 6];
    xlabel('Time (s)', 'FontName', fontName, 'FontSize', fontSize);
    ylabel('Vehicle vertical acceleration (m/s^2)', 'FontName', fontName, 'FontSize', fontSize);
    box on; grid off;
    set(gca, 'FontName', fontName, 'FontSize', fontSize);

    % --- (b) Bridge midspan acceleration (filled marker) ---
    nexttile(2);
    hold on;
    l1 = plot(t_dyn, dynResult.bridge.a, '-s', ...
        'MarkerIndices', idxDyn, 'MarkerFaceColor', CQUBlue, ...
        'Color', CQUBlue, 'LineWidth', lineW);
    lg = legend(l1, legendStr, ...
        'Interpreter', 'latex', 'Location', 'north', ...
        'NumColumns', 1, 'Box', 'off');
    lg.ItemTokenSize = [10, 6];
    xlabel('Time (s)', 'FontName', fontName, 'FontSize', fontSize);
    ylabel('Bridge midspan acceleration (m/s^2)', 'FontName', fontName, 'FontSize', fontSize);
    box on; grid off;
    set(gca, 'FontName', fontName, 'FontSize', fontSize);

    % ---- Figure size ----
    set(gcf, 'Units', 'centimeters', 'Position', [2 2 figW figH]);
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', fontSize);

    %% ================================================================
    %  Figure 4: Bridge coupled frequency (left) & damping ratio (right)
    %% ================================================================
    figure('Name', 'Bridge Modal Variations', 'NumberTitle', 'off');
    tiledlayout(1, 2, 'Padding', 'tight', 'TileSpacing', 'tight');

    % --- (a) Bridge coupled frequency (filled marker) ---
    nexttile(1);
    hold on;
    l1 = plot(eigenResult.t, eigenResult.freq_b, '-s', ...
        'MarkerIndices', idxEig, 'MarkerFaceColor', CQUBlue, ...
        'Color', CQUBlue, 'LineWidth', lineW);
    yline(fb_ref, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
    lg = legend(l1, legendStr, ...
        'Interpreter', 'latex', 'Location', 'north', ...
        'NumColumns', 1, 'Box', 'off');
    lg.ItemTokenSize = [10, 6];
    xlabel('Time (s)', 'FontName', fontName, 'FontSize', fontSize);
    ylabel('Frequency (Hz)', 'FontName', fontName, 'FontSize', fontSize);
    box on; grid off;
    set(gca, 'FontName', fontName, 'FontSize', fontSize);

    % --- (b) Bridge coupled damping ratio (filled marker) ---
    nexttile(2);
    hold on;
    l1 = plot(eigenResult.t, eigenResult.xi_b, '-s', ...
        'MarkerIndices', idxEig, 'MarkerFaceColor', CQUBlue, ...
        'Color', CQUBlue, 'LineWidth', lineW);
    yline(xib_ref, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
    lg = legend(l1, legendStr, ...
        'Interpreter', 'latex', 'Location', 'north', ...
        'NumColumns', 1, 'Box', 'off');
    lg.ItemTokenSize = [10, 6];
    xlabel('Time (s)', 'FontName', fontName, 'FontSize', fontSize);
    ylabel('Damping ratio', 'FontName', fontName, 'FontSize', fontSize);
    box on; grid off;
    set(gca, 'FontName', fontName, 'FontSize', fontSize);

    % ---- Figure size ----
    set(gcf, 'Units', 'centimeters', 'Position', [2 2 figW figH]);
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', fontSize);

end