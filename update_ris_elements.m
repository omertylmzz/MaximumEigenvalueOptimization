function Phi = update_ris_elements(ris_parameters,idx)
    Phi = diag(ris_parameters.beta.*(cos(ris_parameters.phi_rad(:,idx)) + 1i*sin(ris_parameters.phi_rad(:,idx))));
end

