function [bit_stream,reshape_bit_stream,mod_stream,symbol_stream] = stream_generation(sim_parameters)
    bit_stream = randi([0, 1], 1, sim_parameters.nof_bits);
    reshape_bit_stream = reshape(bit_stream, sim_parameters.bps, []).';
    mod_stream = qammod(bi2de(reshape(bit_stream, sim_parameters.bps, []).', 'left-msb'), sim_parameters.mod_order, 'UnitAveragePower', true);
    % mod_stream = mod_stream / sqrt(max(abs(mod_stream).^2));
    symbol_stream = qamdemod(mod_stream, sim_parameters.mod_order);
end

