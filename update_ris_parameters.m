function ris_parameters = update_ris_parameters(ris_parameters,idx)
    ris_parameters.beta = (1 - ris_parameters.beta_min) .* ...
                          (abs(sin(ris_parameters.phi_rad(:,idx) - ris_parameters.c) / 2) .^ ris_parameters.k) + ...
                          ris_parameters.beta_min;
    ris_parameters.A = max(real(ris_parameters.beta));
end

