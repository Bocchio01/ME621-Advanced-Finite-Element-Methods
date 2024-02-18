% ME621 - Advanced Finite Element Methods
% Assignment 2
% Tommaso Bocchietti
% A.Y. 2023/24 - W24

clc
clear variables
close all

%% Problem datas

% Check unit of measure!
rod = struct( ...
    'L0', 200, ...
    'A0', 300, ...
    'E', 70e+9, ...
    'rho0', 2700);

fem = struct( ...
    'N_elements', 5, ...
    'N_nodes_element', 3, ...
    'N_Gauss_points', 3, ...
    'stress_model', "FPK", ...
    'T_final', 0.01, ...
    'F_current', [0 0], ...
    'F_initial', [0 1000]', ...
    'F_final', [1000 1000]');

assert(fem.N_nodes_element >= 2, "The number of nodes per element must be at least 2")


%% Algorithm setup

N_nodes = (fem.N_nodes_element - 1) * fem.N_elements + 1;

step = 0;
time = 0;

u = zeros(N_nodes, 1);
sigma = zeros(N_nodes, 1);

element = struct( ...
    'N0', cell(1), ...
    'B0', cell(1), ...
    'f_int', NaN, ...
    'f_ext', NaN, ...
    'M', NaN);

shape_coefficients = compute_shape_coefficients_1D(fem.N_nodes_element);

N0_handlers = cell(1, fem.N_nodes_element);
B0_handlers = cell(1, fem.N_nodes_element);

for node_element_idx = 1:fem.N_nodes_element
    N0_handlers{node_element_idx} = @(Xi) polyval(shape_coefficients(node_element_idx, :), Xi);
    B0_handlers{node_element_idx} = @(Xi) polyval(polyder(shape_coefficients(node_element_idx, :)), Xi) * 2 / rod.L0;
end

element.N0 = @(Xi) cell2mat(cellfun(@(fun) fun(Xi), N0_handlers, 'UniformOutput', false));
element.B0 = @(Xi) cell2mat(cellfun(@(fun) fun(Xi), B0_handlers, 'UniformOutput', false));

element.f_int = @(rod, fem, u, B0) compute_f_int(rod, fem, u, B0);
element.f_ext = @(rod, fem, u, time) compute_f_ext(rod, fem, time, N0);
element.M = compute_m(rod, fem, element.N0);

clear N0 B0


% Diagonal M

%% Main loop

L = @(e) compute_gather_matrix(e, fem.N_nodes_element, N_nodes);
sum(cellfun(@(e) L(e)' * element.M * L(e), num2cell(1:fem.N_elements), 'UniformOutput', false), 'all')

globalMatrix = @(A) sum(cellfun(@(e) L(e)' * A * L(e), num2cell(1:fem.N_elements), 'UniformOutput', false));

solver = struct( ...
    'f', [], ...
    'M',  diag(sum(globalMatrix(element.M), 2)), ...
    'v', [], ...
    'u', [], ...
    'a', [], ...
    'dt', +inf);

while (time < fem.T_final)
    
    solver.f(step) = 0;
    solver.dt = +inf;
    
    for element_idx = 1:fem.N_elements
        
        u_element = L(element_idx) * u;
        v_element = L(element_idx) * solver.v; % TODO: how to compute v?
        
        f_int_element = 0;
        
        if (step ~= 0)
            f_int_element = element.f_int(rod, fem, u_element, element.B0);
            f_ext_element = element.f_ext(rod, fem, u_element, time);
        end
        
        f_element_kin = f_ext_element - f_int_element;
        
        dt_critical = rod.L0 / sqrt(rod.E / rod.rho0);
        solver.dt = min(solver.dt, dt_critical);
        
        solver.f(step) = solver.f(step) + L(element_idx)' * f_element_kin;
        
    end
    
    solver.dt = 0.875 * solver.dt;
    
    solver.a(step) = solver.M \ solver.f(step);
    
    if (step == 0)
        solver.v(step) = 0 + solver.a(step) * solver.dt;
    else
        solver.v(step) = solver.v(step) + solver.a(step) * solver.dt;
        
    end
    

    solver.u(step+1) = solver.u(step) + solver.v(step) * solver.dt;

    time = time + dt;
    step = step + 1;
    
end




%% Plots

% Xi_vector = -1:0.1:1;
%
%
% nexttile
% hold on
% grid on
%
% for N_idx = 1:length(element.N0)
%     plot(Xi_vector, element.N0(Xi_vector), "LineWidth", 2)
% end
%
% % legend_entries = arrayfun(@(N_idx) ['N', num2str(N_idx)], 1:length(shape_coefficients), 'UniformOutput', false);
% % legend({legend_entries}, 'Location','best')
% title(['Shape functions for a ' num2str(fem.N_nodes_element) '-nodes element'])
%
% nexttile
% hold on
% grid on
%
% for B_idx = 1:length(element.B0)
%     plot(Xi_vector, element.B0{B_idx}(Xi_vector), "LineWidth", 2)
% end
%
% % legend_entries = arrayfun(@(N_idx) ['N', num2str(N_idx)], 1:length(shape_coefficients), 'UniformOutput', false);
% % legend({legend_entries}, 'Location','best')
% title(['Compatibility functions for a ' num2str(fem.N_nodes_element) '-nodes element'])