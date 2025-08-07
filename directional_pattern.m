function F = directional_pattern(alpha, beta)
    % alpha: elevation (0 to pi)
    % beta: azimuth (0 to 2*pi)

    if alpha >= 0 && alpha <= pi/2 && beta >= 0 && beta <= 2*pi
        F = cos(alpha)^3;
    else
        F = 0;
    end
end