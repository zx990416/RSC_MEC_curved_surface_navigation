function MVL = get_mvl(r)
    % Input: r is the firing rate vector for equally spaced directional bins (N x 1)
    
    r = r(:);
    r(r < 0) = 0; 
    
    N = length(r);
    sum_r = sum(r);
    
    if sum_r == 0
        MVL = 0;
        return;
    end

    % Construct equally spaced angle vector theta
    theta = linspace(0, 2*pi, N + 1)';
    theta(end) = []; 
    
    % Compute formula: abs( sum( r_n * e^(-i * theta_n) ) / sum( r_n ) )
    MVL = abs(sum(r .* exp(-1i * theta)) / sum_r);
end