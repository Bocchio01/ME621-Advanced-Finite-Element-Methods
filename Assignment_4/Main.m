% ME621 - Advanced Finite Element Methods
% Assignment 4
% Tommaso Bocchietti
% A.Y. 2023/24 - W24

% vgt_ prefix indicate that the variable is written in Voigt notation
% mat_ prefix indicate that the variable is written in Matrix notation

clc
clear variables
close all


%% Problem datas

data = struct( ...
    'Lx', 1, ...
    'Ly', 1, ...
    'rho', 2700, ...
    'E', 70e+9, ...
    'nu', 0.3, ...
    'vx_top', 0.001);

fem = struct( ...
    'N_nodes_element', 4, ...
    'N_DOF_node', 2, ...
    'N_DOF', NaN, ...
    'Nx', 1, ...
    'Ny', 1, ...
    'constitutive_model', "Cauchy");

output = struct( ...
    'plots', true, ...
    'movie', false);

%% Mesh generation

mesh = struct();
[mesh.coordinates, mesh.connectivity] = generate_mesh(data, fem);
fem.N_DOF = numel(mesh.coordinates);


%% Boundary condition setup

idx_nodes_bottom = reshape((2 * find(mesh.coordinates(:, 2) == 0) + [-1, 0])', 1, []);
idx_nodes_top = reshape((2 * find(mesh.coordinates(:, 2) == data.Ly) + [-1, 0])', 1, []);

initial_coordinates= reshape(mesh.coordinates.', 1, [])';
boundary = struct( ...
    'displacement', struct( ...
    'index', [idx_nodes_bottom idx_nodes_top(2:2:end)], ...
    'value', initial_coordinates([idx_nodes_bottom idx_nodes_top(2:2:end)]) ...
    ), ...
    'velocities', struct( ...
    'index', [idx_nodes_bottom idx_nodes_top(2:2:end) idx_nodes_top(1:2:end)], ...
    'value', [zeros(1, length([idx_nodes_bottom idx_nodes_top(2:2:end)])) data.vx_top * ones(1, length(idx_nodes_top(1:2:end)))] ...
    ));

clear idx_nodes_bottom idx_nodes_top initial_coordinates


%% Funcs structure definitions

funcs = struct( ...
    'integration', struct(), ...
    'mat', struct(), ...
    'vgt', struct(), ...
    'plastic', struct());

% funcs.integration
funcs.integration.weights = [5/9, 8/9, 5/9];
funcs.integration.roots   = [-sqrt(3/5), 0, sqrt(3/5)];
% funcs.integration.weights = [2];
% funcs.integration.roots   = [0];

% funcs.mat
funcs.mat.N      = @(Xi, Eta)    1/4 * [(1-Xi)*(1-Eta), (1+Xi)*(1-Eta), (1+Xi)*(1+Eta), (1-Xi)*(1+Eta)];
funcs.mat.dN     = @(Xi, Eta)    1/4 * [-1+Eta, 1-Eta, 1+Eta, -1-Eta; -1+Xi, -1-Xi, 1+Xi, 1-Xi];
funcs.mat.F      = @(Xi, Eta, u) u * funcs.mat.dN(Xi, Eta)';
funcs.mat.J      = @(Xi, Eta, u) det(funcs.mat.F(Xi, Eta, u));
funcs.mat.B      = @(Xi, Eta, u) inv(funcs.mat.F(Xi, Eta, u))' * funcs.mat.dN(Xi, Eta);
funcs.mat.stress = @(vgt_stress) vgt_stress_to_mat(vgt_stress);
funcs.mat.strain = @(vgt_strain) vgt_strain_to_mat(vgt_strain);

% funcs.vgt
funcs.vgt.N      = @(N)          mat_N_to_vgt(N, fem.N_DOF_node);
funcs.vgt.B      = @(B)          mat_B_to_vgt(B, fem.N_nodes_element, fem.N_DOF_node);
funcs.vgt.stress = @(mat_stress) mat_stress_to_vgt(mat_stress);
funcs.vgt.strain = @(mat_strain) mat_strain_to_vgt(mat_strain);
funcs.vgt.C      = @(E, nu)      E / (1-nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2]; % Plane stress model

% funcs.plastic
funcs.plastic.Sy           = @(S)               sqrt(3/2 * sum(S .* S, 'all'));
funcs.plastic.hardening    = @(strain_p_bar)    200e6 + 325e6 * (strain_p_bar)^0.125;
funcs.plastic.strain_p_bar = @(mat_strain_p)    sqrt(2/3 * sum(mat_strain_p .* mat_strain_p, 'all'));
funcs.plastic.mat_strain_p = @(S, strain_p_bar) 3/2 * strain_p_bar / funcs.plastic.Sy(S) * S;
funcs.plastic.mat_n        = @(S)               3/2 * S / funcs.plastic.Sy(S);
funcs.plastic.phi          = @(S, strain_p_bar) funcs.plastic.Sy(S) - funcs.plastic.hardening(strain_p_bar);


%% States initialization

states = struct( ...
    'dt', [5e-1], ...
    'global', struct(), ...
    'element', struct(), ...
    'plastic', struct(), ...
    'targets', struct());

% states.global
% Here we just initialize the vectors of position and velocities
% Boundary condition are not applied yet
states.global.u(:, 1)      = reshape(mesh.coordinates.', 1, [])';
states.global.v_half(:, 1) = zeros(fem.N_DOF, 1);

states.global.M = zeros(fem.N_DOF, fem.N_DOF);
M = zeros(fem.N_DOF_node * fem.N_nodes_element);

% Here element are the same in initial configuration and we can compute one single element mass matrix knowing
% that it will be the same for all the elements
for i = 1:length(funcs.integration.roots)
    for j = 1:length(funcs.integration.roots)
        
        % Variables for the loop
        idx_nodes = reshape((2 * mesh.connectivity(1, :).' + [-1, 0])', 1, []);
        u = reshape(states.global.u(idx_nodes, 1), [fem.N_DOF_node, fem.N_nodes_element]);
        
        Xi = funcs.integration.roots(i);
        Eta = funcs.integration.roots(j);
        weight_Xi = funcs.integration.weights(i);
        weight_Eta = funcs.integration.weights(j);
        
        N = funcs.mat.N(Xi, Eta);
        J = funcs.mat.J(Xi, Eta, u);
        
        % Loop body
        M = M + weight_Xi * weight_Eta * data.rho * funcs.vgt.N(N)' * funcs.vgt.N(N) * J;
        
    end
end

for element_idx = 1:size(mesh.connectivity, 1)
    
    idx_nodes = reshape((2 * mesh.connectivity(element_idx, :).' + [-1, 0])', 1, []);
    states.global.M(idx_nodes, idx_nodes) = states.global.M(idx_nodes, idx_nodes) + M;
    
end

states.global.M = diag(sum(states.global.M, 2));

% states.element
% Other fields of this struct will be populated directly inside the main loop cycle
states.element.stress = zeros( ...
    size(mesh.connectivity, 1), ...
    length(funcs.integration.roots), ...
    length(funcs.integration.roots), ...
    3);

% states.plastic
states.plastic.strain_p_bar = zeros(size(mesh.connectivity, 1), 1);

% states.targets
% Any variable that require to be stored for each iteration of the main
% lood will be saved in this struct
states.targets.sigma12 = [0];
states.targets.gamma = [0];


clear i j element_idx
clear idx_nodes u Xi Eta weight_Xi weight_Eta N J M


%% Main loop

n = 1;
t = 0;
C = funcs.vgt.C(data.E, data.nu);
dh = 10e-1;

tic
while (states.targets.gamma(n) < 1)
    
    states.global.f(:, n) = zeros(fem.N_DOF, 1);
    % dt_critical = +inf;
    
    % 05) Stress update algorithm
    for element_idx = 1:size(mesh.connectivity, 1)
        
        idx_nodes = reshape((2 * mesh.connectivity(element_idx, :).' + [-1, 0])', 1, []);
        
        states.element.vgt.f_int = zeros(fem.N_DOF_node * fem.N_nodes_element, 1);
        states.element.vgt.f_ext = zeros(fem.N_DOF_node * fem.N_nodes_element, 1);
        
        states.element.vgt.u = states.global.u(idx_nodes, max(n-1, 1));
        states.element.vgt.v_half = states.global.v_half(idx_nodes, max(n-1, 1));
        
        states.element.mat.u = reshape(states.element.vgt.u, [fem.N_DOF_node, fem.N_nodes_element]);
        states.element.mat.v_half = reshape(states.element.vgt.v_half, [fem.N_DOF_node, fem.N_nodes_element]);
        
        if (n > 1)
            for i = 1:length(funcs.integration.roots)
                for j = 1:length(funcs.integration.roots)
                    
                    % Variables for the loop
                    Xi = funcs.integration.roots(i);
                    Eta = funcs.integration.roots(j);
                    weight_Xi = funcs.integration.weights(i);
                    weight_Eta = funcs.integration.weights(j);
                    
                    N = funcs.mat.N(Xi, Eta);
                    B = funcs.mat.B(Xi, Eta, states.element.mat.u);
                    J = funcs.mat.J(Xi, Eta, states.element.mat.u);
                    
                    vgt_stress = squeeze(states.element.stress(element_idx, i, j, :));
                    mat_stress = funcs.mat.stress(vgt_stress);
                    

                    % Elastic predictor
                    switch fem.constitutive_model
                        case "Cauchy"
                            gamma_dot = (states.targets.gamma(n) - states.targets.gamma(n-1)) / states.dt(n-1);
                            
                            F_dot = [1 gamma_dot; 0 1];

                            vgt_strain_rate = funcs.vgt.strain(1/2 * (F_dot + F_dot' - 2*eye(size(F_dot))));
                            
                            vgt_stress_rate = C * vgt_strain_rate;
                            
                        case "Trusdell"
                            mat_L = states.element.mat.v_half * B';
                            vgt_D = funcs.vgt.B(B) * states.element.vgt.v_half;
                            
                            vgt_stress_rate = C * vgt_D + funcs.vgt.stress(mat_L * mat_stress + mat_stress * mat_L' - trace(mat_L) * mat_stress);
                            
                        otherwise
                            error("Constitutive model not implemented yet")
                    end
                    
                    vgt_stress = vgt_stress + vgt_stress_rate * states.dt(n-1);


                    % Plastic corrector
                    mat_stress_dev = compute_stress_deviatoric(vgt_stress);
                    strain_p_bar = states.plastic.strain_p_bar(element_idx);
                    
                    phi = funcs.plastic.phi(mat_stress_dev, strain_p_bar);

                    while(phi > 0)
                        
                        G = data.E / (2 * (1 + data.nu));
                        H = (funcs.plastic.hardening(strain_p_bar + dh) - funcs.plastic.hardening(strain_p_bar)) / dh;
                        % H = 325e6 * 0.125 * (strain_p_bar)^(0.125-1);
                        % error("Having H->+Inf means to have lambda->0 wich create and infite loop (?)")
                        
                        mat_n = funcs.plastic.mat_n(mat_stress_dev);
                        lambda = phi / (3*G + H);
                        
                        vgt_stress = vgt_stress - lambda * C * funcs.vgt.strain(mat_n);
                        
                        strain_p_bar = strain_p_bar + lambda;
                        mat_stress_dev = compute_stress_deviatoric(vgt_stress);
                        
                        phi = funcs.plastic.phi(mat_stress_dev, strain_p_bar);
                        
                    end

                    % Updates elemental forces, stresses and plastic strain 
                    states.element.vgt.f_int = states.element.vgt.f_int + weight_Xi * weight_Eta * funcs.vgt.B(B)' * vgt_stress * J;
                    states.element.vgt.f_ext = states.element.vgt.f_ext + weight_Xi * weight_Eta * funcs.vgt.N(N)' * data.rho * [0; -9.81] * J; % Negligible with respect to the internal forces
                    
                    states.element.stress(element_idx, i, j, :) = vgt_stress;
                    states.plastic.strain_p_bar(element_idx) = strain_p_bar;
                    
                end
            end
        end
        
        states.element.f = states.element.vgt.f_ext - states.element.vgt.f_int;
        states.global.f(idx_nodes, n) = states.global.f(idx_nodes, n) + states.element.f;
        
        % dt_critical = min(compute_min_distance(states.element.mat.u) / sqrt(data.E / data.rho), dt_critical);
        
    end
    
    states.dt(n) = states.dt(1);
    % states.dt(n) = 0.89 * dt_critical;
    
    % 06) Compute acceleration
    states.global.a(:, n) = states.global.M \ states.global.f(:, n);
    
    % 07) Nodal velocities at half-step
    if (n == 1)
        states.global.v_half(:, n) = 0 + 1/2 * states.global.a(:, n) * states.dt(n);
    else
        states.global.v_half(:, n) = states.global.v_half(:, n-1) + states.global.a(:, n) * states.dt(n);
    end
    
    % 08) Enforce velocities BCs
    states.global.v_half(boundary.velocities.index, n) = boundary.velocities.value;
    % if(any(states.targets.gamma > 0.03) && all(states.targets.sigma12 >= 0))
    %     states.global.v_half(boundary.velocities.index, n) = -boundary.velocities.value;
    % end 
    % 
    % if(any(states.targets.gamma > 0.05) && all(states.targets.gamma >= -0.02))
    %     states.global.v_half(boundary.velocities.index, n) = -boundary.velocities.value;
    % end 

    % 08) Check energy balance
    % funcs.w_int(:, n) = funcs.w_int(:, max(n-1, 1)) + (funcs.u(:, n) - funcs.u(:, max(n-1, 1)))' * ((funcs.f(:, n) - funcs.f(:, max(n-1, 1)))/2);
    % funcs.w_ext(:, n) = funcs.w_ext(:, max(n-1, 1)) + (funcs.u(:, n) - funcs.u(:, max(n-1, 1)))' * ((funcs.f(:, n) - funcs.f(:, max(n-1, 1)))/2);
    % funcs.w_kin(:, n) = 1/2 * funcs.v_half(:, n)' * funcs.M * funcs.v_half(:, n);
    % funcs.w(n) = 1e-3 * (1e-3 * funcs.w_kin(:, n) + funcs.w_int(:, n) - funcs.w_ext(:, n));
    % assert(abs(funcs.w(:, n) - funcs.w(:, max(n-1, 1))) < eps, "Energy Balance not verified")
    
    % 09) Update nodal displacements
    if (n == 1)
        states.global.u(:, n) = states.global.u(:, n) + states.global.v_half(:, n) * states.dt(n);
    else
        states.global.u(:, n) = states.global.u(:, n-1) + states.global.v_half(:, n) * states.dt(n);
    end
    
    % 10) Enforce displacements BCs
    states.global.u(boundary.displacement.index, n) = boundary.displacement.value;
    
    % 11) Update time_step n and time t
    states.targets.sigma12(n + 1) = states.element.stress(1, 1, 1, 3);
    states.targets.gamma(n + 1) = atan(states.global.u(3, n) / states.global.u(4, n));
    
    t = t + states.dt(n);
    n = n + 1;
    
    if (output.movie && mod(n, 10) == 0)
        set(0,'DefaultFigureVisible','off');
        if (exist('plot_figures', 'var') == 0)
            plot_figures(1) = getframe(plot_states(mesh, states, t));
        end
        plot_figures(end+1) = getframe(plot_states(mesh, states, t));
    end
    
end
toc

clear i j element_idx idx_nodes
clear n t C G H lambda mat_n mat_stress_dev strain_p_bar vgt_strain_rate phi dh dt_critical
clear Xi Eta weight_Xi weight_Eta N B J vgt_stress mat_stress mat_L vgt_D vgt_stress_rate vgt_stress


%% Movie

if (output.movie)
    
    video_object = VideoWriter('FEM.avi', 'Motion JPEG AVI');
    video_object.Quality = 80;
    video_object.FrameRate = 10;
    
    open(video_object);
    writeVideo(video_object, plot_figures);
    close(video_object);
    
end


%% Plots

set(0,'DefaultFigureVisible','on');
if output.plots
    
    figure_results_analysis = plot_states(mesh, states);
    
end

clear figure_results_analysis
