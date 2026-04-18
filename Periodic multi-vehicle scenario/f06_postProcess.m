%% f06_postProcess.m
% =========================================================================
% Description : Post-processing and visualization for the periodic
%               multi-vehicle VBI system.
%               Generates two figures:
%               Fig.1 - Bridge mid-span displacement, velocity, acceleration
%               Fig.2 - Coupled bridge frequency and damping ratio vs time
% Author      : Shi LuoKe (史罗克)
% Date        : 2026-04-18
% -------------------------------------------------------------------------
% Input  : eigenResult - Coupled eigenvalue results
%          dynResult   - Dynamic response results
%          theorFreq   - Theoretical uncoupled frequencies
%          sysParams   - System parameter structure
% =========================================================================

function f06_postProcess(eigenResult, dynResult, theorFreq, sysParams)

    %% Unpack
    dt        = sysParams.solver.dt;
    gamma_val = sysParams.solver.gamma;
    Num_steps = length(dynResult.bridge.d);
    t_dyn     = (0:Num_steps-1) * dt;

    % Theoretical reference values
    fb_ref  = theorFreq.bridge.f1;
    xib_ref = sysParams.bridge.xib;

    % ---- Common style settings ----
    CQUBlue   = '#02529f';
    fontName  = 'Times New Roman';
    fontSize  = 11;
    lineW     = 1;
    figW      = 14;
    figH      = 6;

    % Legend string
    legendStr = sprintf('$\\gamma = %g$', gamma_val);

    % Marker indices
    numMarkers = 8;
    idxDyn = round(linspace(1, Num_steps, numMarkers));
    idxEig = round(linspace(1, length(eigenResult.t), numMarkers));



%% ================================================================
    %  Figure: Bridge coupled frequency
    %% ================================================================
    figure('Name', 'Bridge Modal Variations', 'NumberTitle', 'off');

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
    set(gcf, 'Units', 'centimeters', 'Position', [2 2 7 6]);
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', fontSize);
end