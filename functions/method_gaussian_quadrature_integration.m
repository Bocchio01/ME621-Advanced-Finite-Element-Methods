% roots = vpasolve(legendreP(7,x) == 0)

function result = method_gaussian_quadrature_integration(integrand_func, fem)

switch fem.N_Gauss_points
    case 1
        weights = [2];
        roots = [0];
    case 2
        weights = [1, 1];
        roots = [-1/sqrt(3), 1/sqrt(3)];
    case 3
        weights = [5/9, 8/9, 5/9];
        roots = [-sqrt(3/5), 0, sqrt(3/5)];
    otherwise
        error('Order of Gaussian quadrature not implemented (too high)')
end

result = 0;
for point_idx = 1:length(roots)
    
    result = result + integrand_func(roots(point_idx)) * weights(point_idx);
    
end

end