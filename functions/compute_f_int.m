function f_int = compute_f_int(rod, fem, u, B0)


switch fem.stress_model
    case "FPK"
        integrand_func = @(Xi) B0(Xi)' * rod.E * B0(Xi) * u * rod.A0 * rod.L0 / 2;
        
    otherwise
        error("Stress model not implemented")
        
end

f_int = method_gaussian_quadrature_integration(integrand_func, fem);

end