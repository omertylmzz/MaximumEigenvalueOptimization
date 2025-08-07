function [PL,noise] = pl_noise_generation(sim_parameters,ris_parameters)
    switch sim_parameters.pl_model
        case "ETSI Near Field"
            PL = ris_parameters.Gt*ris_parameters.Gr*((sim_parameters.lambda^2)/(16*pi^2)) ...
                    * ((ris_parameters.A/(ris_parameters.dt+ris_parameters.dr))^2);
        case "ETSI Far Field"
            PL = (64*pi^3*(ris_parameters.dt*ris_parameters.dr)^2) ...
                    / (ris_parameters.Gt*ris_parameters.Gr*ris_parameters.Gs*(sim_parameters.nof_ris_elements)^2 ...
                       *ris_parameters.dx*ris_parameters.dy*sim_parameters.lambda^2 ...
                       *ris_parameters.Ftx*ris_parameters.Frx*ris_parameters.A^2);
    end
    noise = (sim_parameters.N0) * (randn(sim_parameters.nof_rx_antennas,1) + 1i*randn(sim_parameters.nof_rx_antennas,1));
end

