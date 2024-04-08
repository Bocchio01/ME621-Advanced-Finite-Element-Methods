function [delta_strain_plastic, strain_new] = method_radial_return(E, nu, H, Sy0, strain, strain_plastic, Deps)

mu = E / (2 * (1 + nu));
kappa = E / (3 * (1 - 2*nu));
vl=[1 1 1 0]';

trDeps = Deps(1) + Deps(2);
De =  [Deps(1:2); 0; Deps(3)/2] - 1/3 * trDeps * vl;

strain_elas = strain + kappa * trDeps * vl + 2 * mu * De;
S = strain_elas - 1/3 * sum(strain_elas(1:3)) * vl;

Sy_vonMises = sqrt(1.5 * (S(1)^2 + S(2)^2 + S(3)^2 + 2*S(4)^2) );

phi = Sy_vonMises - (Sy0 + H * strain_plastic);

if (phi > 0)
    
    delta_strain_plastic = phi / (3*mu+H);
    sigeq_new = Sy_vonMises - 3*mu*delta_strain_plastic;
    Depsp = 3/2 * delta_strain_plastic * S/Sy_vonMises;
    strain_new = strain_elas - 2*mu*Depsp;
    
else
    
    delta_strain_plastic = 0;
    strain_new = strain_elas;
    
end

end