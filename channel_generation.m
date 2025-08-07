function [H,G] = channel_generation(sim_parameters)
    if ~sim_parameters.dp_channel_model
        H = sqrt(.5)*(randn(sim_parameters.nof_ris_elements, sim_parameters.nof_tx_antennas,sim_parameters.nof_symbols) ...
                      + 1i*randn(sim_parameters.nof_ris_elements, sim_parameters.nof_tx_antennas,sim_parameters.nof_symbols));
        G = sqrt(.5)*(randn(sim_parameters.nof_rx_antennas, sim_parameters.nof_ris_elements, sim_parameters.nof_symbols) ...
                      + 1i*randn(sim_parameters.nof_rx_antennas, sim_parameters.nof_ris_elements, sim_parameters.nof_symbols));
    else
        dp_channel_params = sim_parameters.dp_channel_params;
        psi = pi*dp_channel_params.a_psi/2;
        phi = pi*dp_channel_params.a_phi/2;
        vx = exp(1j*dp_channel_params.dlx*2*pi*(0:(sim_parameters.nof_ant_ports_N1)-1)'*(sin(psi).*cos(phi)));
        vy = exp(1j*dp_channel_params.dly*2*pi*(0:(sim_parameters.nof_ant_ports_N2)-1)'*(sin(psi).*sin(phi)));
        At  = kron(vy,vx);   

        doa = pi*dp_channel_params.a_doa*2/3;
        Ar = exp(1j*dp_channel_params.dlr*[0:sim_parameters.nof_ant_ports_Nr-1]'*2*pi*sin(doa));

        alpha = sqrt(size(At,1)*sim_parameters.nof_ant_ports_Nr) * ...
                (sqrt(.5)*(randn(sim_parameters.nof_ris_elements, sim_parameters.nof_tx_antennas,sim_parameters.nof_symbols) ...
                      + 1i*randn(sim_parameters.nof_ris_elements, sim_parameters.nof_tx_antennas,sim_parameters.nof_symbols)));
        
        
        for i1 = 1:2
            for i2 = 1:2
                beta = alpha.*exp(1j*pi*randn(nPaths,1));
                beta(1) = beta(1)*sqrt(kappa/(kappa+1));
                beta(2:end) = beta(2:end)*sqrt(1/(kappa+1));
                B = [B,beta];
                H((i1-1)*Mr+1:i1*Mr, (i2-1)*NN+1:i2*NN) = sqrt(kappa/(kappa+1))*Ar * diag(beta) * At';
            end
        end

    end
end

