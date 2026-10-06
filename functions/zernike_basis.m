function Z = zernike_basis(n, m, r, theta)

    m_abs = abs(m);
    if mod(n - m_abs, 2) ~= 0
        error('Invalid Zernike degree combination: n - |m| must be even.');
    end

    R = zeros(size(r));
    for s = 0:((n - m_abs) / 2)
        num = ((-1)^s) * factorial(n - s);
        den = factorial(s) * factorial((n + m_abs)/2 - s) * factorial((n - m_abs)/2 - s);
        R = R + (num / den) * (r.^(n - 2 * s));
    end

    if m >= 0
        Z = R .* cos(m * theta);
    else
        Z = R .* sin(m_abs * theta);
    end
end
