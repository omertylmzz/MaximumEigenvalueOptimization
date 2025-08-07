clc; clear; close all;
disp("##############################################################################################")
%% SVD-Based
[sim_parameters,ris_parameters] = load_parameters("SVD-Based");
[bit_stream,reshape_bit_stream,mod_stream,symbol_stream] = stream_generation(sim_parameters);
[H,G] = channel_generation(sim_parameters);
%% Main Loop
for snr_idx = 1:numel(sim_parameters.snr_range)
    Es = 10^(sim_parameters.snr_range(snr_idx)/10)*sim_parameters.nof_tx_antennas*sim_parameters.N0^2;
    berror = 0;
    serror = 0;
    power_rx = 0;
    for idx = 1:sim_parameters.nof_symbols  
        % Update RIS Parameters & Configuration
        ris_parameters = update_ris_parameters(ris_parameters,idx);
        [PL,noise] = pl_noise_generation(sim_parameters,ris_parameters);
        Phi = update_ris_elements(ris_parameters,idx);
        % MIMO-RIS Channel & Precoding Vector Generation
        F = G(:,:,idx) * Phi * H(:,:,idx); 
        [~,V,U] = svd(F);
        if sim_parameters.codebook_type == "SVD-Based" 
            W = U(:,1);
        else
            W = selection_precoder(sim_parameters,U);
        end
        % Transceiver Process (Including BER & SER Calculation)
        tx_stream = W*mod_stream(idx);
        rx_stream = F*sqrt(PL)*Es*tx_stream + noise;
        wzf = pinv((F*W)'*F*W)*(F*W)'/(sqrt(PL)*Es);
        equalized_stream = wzf*rx_stream;
        power_rx = power_rx + (abs(wzf * (F * sqrt(PL) * Es * tx_stream))^2 / abs(wzf * noise)^2);
        rx_bit_stream = qamdemod(equalized_stream, sim_parameters.mod_order,'UnitAveragePower',true,OutputType='bit').';
        rx_symbol_stream = qamdemod(equalized_stream, sim_parameters.mod_order);
        berror = berror + sum(reshape_bit_stream(idx,:) ~= rx_bit_stream);
        serror = serror + (symbol_stream(idx) ~= rx_symbol_stream);
    end
    sim_parameters.berror(snr_idx) = (berror / sim_parameters.nof_bits)*100;
    sim_parameters.serror(snr_idx) = (serror / sim_parameters.nof_symbols)*100;
    sim_parameters.power_rx(snr_idx) = 10 * log10(power_rx/sim_parameters.nof_symbols);
    fprintf('Progress: %3.0f%% [%d/%d]   | BER: %3.0f%% | SER: %3.0f%% | Equalized SNR (dB): %3.0f \n', ...
            (snr_idx/numel(sim_parameters.snr_range))*100, ...
            snr_idx, ...
            numel(sim_parameters.snr_range), ...
            sim_parameters.berror(snr_idx), ...
            sim_parameters.serror(snr_idx), ...
            sim_parameters.power_rx(snr_idx));
end

disp("##############################################################################################")
clear;
%% Type1
[sim_parameters,ris_parameters] = load_parameters("Type1");
[bit_stream,reshape_bit_stream,mod_stream,symbol_stream] = stream_generation(sim_parameters);
[H,G] = channel_generation(sim_parameters,ris_parameters);
%% Main Loop
for snr_idx = 1:numel(sim_parameters.snr_range)
    Es = 10^(sim_parameters.snr_range(snr_idx)/10)*sim_parameters.nof_tx_antennas*sim_parameters.N0^2;
    berror = 0;
    serror = 0;
    power_rx = 0;
    for idx = 1:sim_parameters.nof_symbols  
        % Update RIS Parameters & Configuration
        ris_parameters = update_ris_parameters(ris_parameters,idx);
        [PL,noise] = pl_noise_generation(sim_parameters,ris_parameters);
        Phi = update_ris_elements(ris_parameters,idx);
        % MIMO-RIS Channel & Precoding Vector Generation
        F = G(:,:,idx) * Phi * H(:,:,idx); 
        [~,V,U] = svd(F);
        if sim_parameters.codebook_type == "SVD-Based" 
            W = U(:,1);
        else
            W = selection_precoder(sim_parameters,U);
        end
        % Transceiver Process (Including BER & SER Calculation)
        tx_stream = W*mod_stream(idx);
        rx_stream = F*sqrt(PL)*Es*tx_stream  + noise;
        wzf = pinv((F*W)'*F*W)*(F*W)'/(sqrt(PL)*Es);
        equalized_stream = wzf*rx_stream;
        power_rx = power_rx + (abs(wzf * (F * sqrt(PL) * Es * tx_stream))^2 / abs(wzf * noise)^2);
        rx_bit_stream = qamdemod(equalized_stream, sim_parameters.mod_order,'UnitAveragePower',true,OutputType='bit').';
        rx_symbol_stream = qamdemod(equalized_stream, sim_parameters.mod_order);
        berror = berror + sum(reshape_bit_stream(idx,:) ~= rx_bit_stream);
        serror = serror + (symbol_stream(idx) ~= rx_symbol_stream);
    end
    sim_parameters.berror(snr_idx) = (berror / sim_parameters.nof_bits)*100;
    sim_parameters.serror(snr_idx) = (serror / sim_parameters.nof_symbols)*100;
    sim_parameters.power_rx(snr_idx) = 10 * log10(power_rx/sim_parameters.nof_symbols);
    fprintf('Progress: %3.0f%% [%d/%d]  | BER: %3.0f%% | SER: %3.0f%% | Equalized SNR (dB): %3.0f \n', ...
            (snr_idx/numel(sim_parameters.snr_range))*100, ...
            snr_idx, ...
            numel(sim_parameters.snr_range), ...
            sim_parameters.berror(snr_idx), ...
            sim_parameters.serror(snr_idx), ...
            sim_parameters.power_rx(snr_idx));
end



