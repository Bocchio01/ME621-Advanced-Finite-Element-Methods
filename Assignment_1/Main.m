% ME621 - Advanced Finite Element Methods
% Assignment 1
% Tommaso Bocchietti
% A.Y. 2023/24 - W24
% 
% Models obtained using Mathematica. Refer to the attached Report.pdf 
% for further details.

clc
clear variables
close all


%% Problem datas

geometry = struct( ...
    'L', [3 0.5], ...
    'A', 1e-3 * [1 1], ...
    'E', 70e+9 * [1 1]);

algorithm = struct( ...
    'N_steps', 1000, ...
    'U_zero', [0 0]', ...
    'P_zero', [0 0]', ...
    'P_delta', 5 * [1 1]', ...
    'tolerance', 1e-5, ...
    'correctors', ["NR", "mNR"]);

models = {
    @(U, geometry) model_exact(U, geometry), ...
    @(U, geometry) model_taylor_1(U, geometry), ...
    @(U, geometry) model_taylor_2(U, geometry), ...
    @(U, geometry) model_taylor_3(U, geometry), ...
    @(U, geometry) model_approximated_soft(U, geometry),...
    @(U, geometry) model_approximated_hard(U, geometry),...
    };

output = struct( ...
    'plots', true * [1 1 1], ...
    'export',  false);

assert(isequal(size(geometry.L), size(geometry.A), size(geometry.E)), 'Inconsistent input geometry data');


%% Solver structure
% Used to store all the models and results of the script.
% Organized as follow:
% - model_name: name of the model
% - handler:    refers to the anonymous function that returns [Kt, Fi]
%               given the current position and geometry [U, geometry]
% - results:    struct
%       - corrector#1 -> results obtained using corrector#1
%       - corrector#2 -> results obtained using corrector#2
%       - ...
%       - corrector#N -> results obtained using corrector#N
% 
% Each results entry has the following structure:
% - corrector:  name of the corrector used to solve
% - U:          vector/matrix containing displachment for each step (N_steps x N_displachments)
% - iteraction: vector containing number of iterations needed for each step (N_steps x 1)
% - time:       vector containing CPU time used for each step (N_steps x 1)

solvers = repmat( ...
    create_solver_struct(length(algorithm.U_zero), algorithm.N_steps, length(algorithm.correctors)), ...
    length(models), ...
    1);

for i = 1:length(models)

    solvers(i).model_name = strrep(func2str(models{i}), '_', '-');
    solvers(i).handler = models{i};

    for j = 1:length(algorithm.correctors)
        solvers(i).results(j).corrector = algorithm.correctors{j};
    end

end

clear models corrector_name i j


%% Models solution

inverse = @(A) A \ eye(size(A));
residual = @(P, Fi) Fi - P;

for solver_idx = 1:length(solvers)

    for corrector_idx = 1:length(algorithm.correctors)

        U = algorithm.U_zero;
        P = algorithm.P_zero;

        for step = 1:algorithm.N_steps

            P = P + algorithm.P_delta;

            % Start predictor
            [Kt, Fi] = solvers(solver_idx).handler(U, geometry);

            Kt0_inv = inverse(Kt);
            U = method_euler(U, algorithm.P_delta, Kt0_inv);
            % End predictor

            tic;
            iteractions = 0;
            % Start corrector
            while(norm(residual(P, Fi), 1) > algorithm.tolerance)

                iteractions = iteractions + 1;

                [Kt, Fi] = solvers(solver_idx).handler(U, geometry);

                switch corrector_idx
                    case 1
                        Kt_inv = inverse(Kt);
                    case 2
                        Kt_inv = Kt0_inv;
                    otherwise
                        error('Unknown corrector name');
                end

                U = method_newton_raphson(U, residual(P, Fi), Kt_inv);

            end
            % End corrector
            time = toc;

            solvers(solver_idx).results(corrector_idx).U(step, :) = U;
            solvers(solver_idx).results(corrector_idx).iterations(step) = iteractions + 1;
            solvers(solver_idx).results(corrector_idx).time(step) = time;

        end
    end
end

clear model_names model_name corrector_idx corrector_name step
clear P U time iteractions Kt Fi Kt_inv Kt0_inv inverse residual


%% Export csv

% numerical_results_Pxxxx_Pyxxx
if output.export
    
    N_correctors = length(algorithm.correctors);
    N_models = length(solvers);

    model_names = reshape(repmat({solvers.model_name}, 1, 1), 1, [])';

    P_max = algorithm.N_steps * algorithm.P_delta + algorithm.P_zero;
    P_maxs = reshape(repmat(P_max, N_models, 1), 2, [])';
    U_maxs = zeros(N_models, 2);
    
    for solver_idx = 1:N_models
        U_maxs(solver_idx, :) = solvers(solver_idx).results(1).U(end, :);
    end

    % Convert numeric arrays to cell arrays
    model_names_cell = cellstr(model_names);
    P_maxs_cell = num2cell(P_maxs);
    U_maxs_cell = arrayfun(@(x) sprintf('%.2e', x), U_maxs, 'UniformOutput', false);

    % Concatenate cell arrays
    data = [model_names_cell, P_maxs_cell, U_maxs_cell];

    header = ["Model name" "Px" "Py" "Ux" "Uy"];
    writematrix([header; data], 'Assignment_1/latex/files/numerical_results_Px5000000_Py5000000.csv');

end


%% Plots

points = [
    0, 0.5; ...
    3, 0.5; ...
    3, 0];

if output.plots(1)

    figure_deformed_shape = figure('Name', 'Deformed shape analysis', 'NumberTitle', 'off');
    tiledlayout(2, 3);

    nexttile([1 3]);
    scale_factor = 1e1;
    hold on
    grid on

    axis equal;
    % xlim([-0.5, 5]);
    % ylim([-0.1, 1]);

    plot(points(:, 1), points(:, 2), '-k', 'LineWidth', 2);

    for solver_idx = 1:length(solvers)

        B = scale_factor * solvers(solver_idx).results(1).U(end, :) + [3 0.5];
        pointsDeformed = [
            points(1, :);
            B;
            points(3, :)
            ];
        plot(pointsDeformed(:, 1), pointsDeformed(:, 2), '--', 'LineWidth', 1);

    end

    title(sprintf("Final deformed shape @scaleFactor = %d", scale_factor))
    legend(["Unloaded structure"; {solvers.model_name}'], 'Location','best')
    xlabel('x [m]');
    ylabel('y [m]');


    % Trajectory of the deformed shape
    nexttile(4);
    scale_factor = 1e6;
    hold on
    grid on
    axis equal
    % xlim([0, scale_factor * max(solvers(solver_idx).results(1).U(:, 1))]);
    % ylim([0, scale_factor * max(solvers(solver_idx).results(1).U(:, 2))]);

    if (size(algorithm.U_zero, 1) == 2)
        for solver_idx = 1:length(solvers)

            B = solvers(solver_idx).results(1).U;
            plot(scale_factor * [0; B(:, 1)], ...
                scale_factor * [0; B(:, 2)], ...
                '-', 'LineWidth', 1);

        end
    end

    title("Trajectory of B point")
    legend({solvers.model_name}, 'Location','best')
    xlabel('u [\mum]');
    ylabel('v [\mum]');


    % Error of Ux with respect to exact model
    nexttile(5);
    scale_factor = 1e6;
    hold on
    grid on

    for solver_idx = 1:length(solvers)

        plot(algorithm.P_delta(1) * 0:algorithm.N_steps-1, ...
            scale_factor * (solvers(solver_idx).results(1).U(:, 1) - solvers(1).results(1).U(:, 1)), ...
            '-', 'LineWidth', 1);

    end

    title("Error for U_x")
    legend({solvers.model_name}, 'Location','best')
    xlabel('P_x [Pa]');
    ylabel('E_u [\mum]');


    % Error of Uy with respect to exact model
    nexttile(6);
    scale_factor = 1e6;
    hold on
    grid on

    for solver_idx = 1:length(solvers)

        plot(algorithm.P_delta(2) * 0:algorithm.N_steps-1, ...
            scale_factor * (solvers(solver_idx).results(1).U(:, 2) - solvers(1).results(1).U(:, 2)), ...
            '-', 'LineWidth', 1);

    end

    title("Error for U_y")
    legend({solvers.model_name}, 'Location','best')
    xlabel('P_y [Pa]');
    ylabel('E_v [\mum]');

    figs(1) = figure_deformed_shape;

    clear solver_idx scale_factor B figure_deformed_shape pointsDeformed
end

if output.plots(2)

    figure_correctors_comparison = figure('Name', 'Correctors methods comparison', 'NumberTitle', 'off');
    tiledlayout(2, 1);

    % Iteraction comparison
    nexttile(1);
    hold on
    grid on

    % xlim([0.5, length(solvers) + 0.5]);
    % ylim([0, 5]);

    xticks(1:length(solvers));
    xticklabels({solvers.model_name});

    NR = zeros(length(solvers), 1);
    mNR = zeros(length(solvers), 1);

    for solver_idx = 1:length(solvers)

        NR(solver_idx) = mean(solvers(solver_idx).results(1).iterations);
        mNR(solver_idx) = mean(solvers(solver_idx).results(2).iterations);

    end

    plot(NR, '-*', 'LineWidth', 2);
    plot(mNR, '-*', 'LineWidth', 2);

    title('Iteraction comparison')
    legend(algorithm.correctors, 'Location','best')
    xlabel('Models');
    ylabel('#iteraction');


    % Time comparison
    nexttile(2);
    scale_factor = 1e6;
    hold on
    grid on

    % xlim([0.5, length(solvers) + 0.5]);
    % ylim([0, 70]);

    xticks(1:length(solvers));
    xticklabels({solvers.model_name});

    NR = zeros(length(solvers), 1);
    mNR = zeros(length(solvers), 1);

    for solver_idx = 1:length(solvers)

        NR(solver_idx) = mean(solvers(solver_idx).results(1).time);
        mNR(solver_idx) = mean(solvers(solver_idx).results(2).time);

    end

    plot(scale_factor * NR, '-*', 'LineWidth', 2);
    plot(scale_factor * mNR, '-*', 'LineWidth', 2);

    title('Time comparison')
    legend(algorithm.correctors, 'Location','best')
    xlabel('Models');
    ylabel('Computational time [\mus]');

    figs(2) = figure_correctors_comparison;

    clear solver_idx scale_factor NR mNR figure_correctors_comparison

end

if output.plots(3)

    figure_force_displacement = figure('Name', 'Force versus displacement', 'NumberTitle', 'off');
    tiledlayout(1, 2);

    % Force(displacement) along x direction
    nexttile(1);
    scale_factor = 1e6;
    hold on
    grid on

    for solver_idx = 1:length(solvers)

        plot(scale_factor * [0; solvers(solver_idx).results(1).U(:, 1)], ...
            algorithm.P_delta(1) * (0:algorithm.N_steps), ...
            '-', 'LineWidth', 1);

    end

    title('P_x(U_x)')
    legend({solvers.model_name}, 'Location','best')
    xlabel('U_x [\mum]');
    ylabel('P_x [Pa]');

    % Force(displacement) along y direction
    nexttile(2);
    scale_factor = 1e6;
    hold on
    grid on

    for solver_idx = 1:length(solvers)

        plot(scale_factor * [0; solvers(solver_idx).results(1).U(:, 2)], ...
            algorithm.P_delta(2) * (0:algorithm.N_steps), ...
            '-', 'LineWidth', 1);

    end

    title('P_y(U_y)')
    legend({solvers.model_name}, 'Location','best')
    xlabel('U_y [\mum]');
    ylabel('P_y [Pa]');

    figs(3) = figure_force_displacement;

    clear solver_idx scale_factor figure_force_displacement

end


for fig_idx = sort(find(output.plots == true), 'descend')

    figure(figs(fig_idx));
    
    if output.export
        figureName = lower(strrep(get(figs(fig_idx), 'Name'), ' ', '_'));
        % set(figs(fig_idx), 'Position', [100, 100, 800, 600]);
        set(figs(fig_idx), 'PaperPosition', 2 * [0, 0, 15, 10]);
        print(figs(fig_idx), ['Assignment_1\Results\' figureName '.png'], '-dpng', '-r100')
    end
    
end

clear points fig_idx


%% Functions

function solver_struct = create_solver_struct(U_size, N_steps, N_correctors)

result_struct = create_result_struct(U_size, N_steps);

solver_struct = struct( ...
    'model_name', '', ...
    'handler', @handle, ...
    'results', repmat(result_struct, N_correctors, 1) ...
    );

end

function result_struct = create_result_struct(U_size, N_steps)

result_struct = struct( ...
    'corrector', '', ...
    'U', zeros(N_steps, U_size), ...
    'iterations', zeros(N_steps, 1), ...
    'time', zeros(N_steps, 1) ...
    );

end







