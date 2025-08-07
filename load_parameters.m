function [sim_parameters,ris_parameters] = load_parameters(codebook_type)
    % General Parameters
    sim_parameters.nof_ant_ports_XP = 2;
    sim_parameters.nof_ant_ports_N1 = 4;
    sim_parameters.nof_ant_ports_N2 = 2;
    sim_parameters.nof_layers = 1;
    sim_parameters.nof_tx_antennas = sim_parameters.nof_ant_ports_XP*sim_parameters.nof_ant_ports_N1*sim_parameters.nof_ant_ports_N2;
    sim_parameters.nof_ant_ports_Nr = 8;
    sim_parameters.nof_rx_antennas = sim_parameters.nof_ant_ports_Nr*sim_parameters.nof_ant_ports_XP;
    sim_parameters.nof_ris_elements = 100;
    sim_parameters.snr_range = 25:50;
    sim_parameters.berror = zeros(numel(sim_parameters.snr_range),1);
    sim_parameters.serror = zeros(numel(sim_parameters.snr_range),1);
    sim_parameters.power_rx = zeros(numel(sim_parameters.snr_range),1);
    sim_parameters.dp_channel_model = true;
    if sim_parameters.dp_channel_model
        sim_parameters.dp_channel_params.dlx = 0.5;
        sim_parameters.dp_channel_params.dly = 0.5;
        sim_parameters.dp_channel_params.dlr = 0.5;

        sim_parameters.dp_channel_params.a_psi = rand(1);
        sim_parameters.dp_channel_params.a_phi = rand(1)-0.5;
        sim_parameters.dp_channel_params.a_doa = rand(1)-0.5;

        sim_parameters.dp_channel_params.K = 0;
    end
    sim_parameters.N0 = 10^(-1.5);
    sim_parameters.pl_model = "ETSI Near Field";
    sim_parameters.codebook_type = codebook_type; % "SVD-Based" or "Type1"
    sim_parameters.codebook_mode = 1;
    if sim_parameters.codebook_type == "Type1"
        sim_parameters.codebook = generate_codebook(sim_parameters);
    end
    
    % Scheduler Parameters
    sim_parameters.mod_order = 16;
    sim_parameters.bps = log2(sim_parameters.mod_order);
    sim_parameters.nof_bits = 1e6;
    sim_parameters.nof_symbols = sim_parameters.nof_bits / sim_parameters.bps;
    
    % Environment Parameters
    sim_parameters.fc = 3.6e9;
    sim_parameters.c0 = physconst('LightSpeed');
    sim_parameters.lambda = sim_parameters.c0 / sim_parameters.fc;
    
    % RIS Parameters
    ris_parameters.beta_min = 0.8;
    ris_parameters.c = 0.43 * pi;
    ris_parameters.k = 1.6;
    ris_parameters.qlevels = [0, 90, 180, 270]; 
    ris_parameters.phi_rad = deg2rad(ris_parameters.qlevels(randi(length(ris_parameters.qlevels), ...
                                                                  sim_parameters.nof_ris_elements, ...
                                                                  sim_parameters.nof_symbols)));  
    ris_parameters.beta = zeros(sim_parameters.nof_ris_elements,1);
    ris_parameters.A = 0;
    
    ris_parameters.dx = sim_parameters.lambda / 4;
    ris_parameters.dy = sim_parameters.lambda / 4;
    
    dtx = 30; 
    dty = 40;
    drx = dtx; 
    dry = dty;
    ris_parameters.dt = sqrt(dtx^2+dty^2); 
    ris_parameters.dr = sqrt(drx^2+dry^2); 
    
    ris_parameters.Gt = 1; 
    ris_parameters.Gr = 1;
    ris_parameters.Gs = (4*pi*ris_parameters.dx*ris_parameters.dy)/sim_parameters.lambda^2;
    
    theta_t = deg2rad(0); 
    phi_t = deg2rad(150);
    theta_r = deg2rad(0); 
    phi_r = deg2rad(30);
    ris_parameters.Ftx = directional_pattern(theta_t, phi_t);
    ris_parameters.Frx = directional_pattern(theta_r, phi_r);
end

