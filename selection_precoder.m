function W = selection_precoder(sim_parameters, U)
    W = zeros(size(U(:,1)));
    [i11,i12] = selection_i11_i12(sim_parameters.codebook, U(1:floor(end/2),1));
    if sim_parameters.nof_layers > 1
        i13 = selection_i13(sim_parameters.codebook,U(1:floor(end/2),2),i11,i12);
    end
    i2 = selection_i2(sim_parameters.codebook,U(:,1),i11,i12);
    if sim_parameters.nof_layers == 1
        W = sim_parameters.codebook(:,:,i2,i11,i12);
    else
        W = sim_parameters.codebook(:,:,i2,i11,i12,i13).';
    end
end

function [i11,i12] = selection_i11_i12(codebook, wopt)
    corr = 0; 
    max_corr = 0;
    for i = 1:size(codebook,4)
        for j = 1:size(codebook,5)
            corr = abs(codebook(1:size(codebook,1)/2,1,1,i,j)'*wopt);
            if max_corr < corr; max_corr = corr; i11 = i; i12 = j; end
        end
    end
end

function i13 = selection_i13(codebook,wopt,i11,i12)
    corr = 0; 
    max_corr = 0;
    for i = 1:size(codebook,6)
        corr = abs(codebook(1:size(codebook,1)/2,2,1,i11,i12,i)'*wopt);
        if max_corr < corr; max_corr = corr; i13 = i; end
    end
end

function i2 = selection_i2(codebook,wopt,i11,i12)
    corr = zeros(size(codebook,3),1);
    vlm = codebook(1:size(codebook,1)/2,1,1,i11,i12).';
    for i = 1:size(codebook,3)
        temp = [vlm (exp(1i*pi*(i-1)/2) * vlm)].';
        corr(i) = abs(temp'*wopt);
    end
    [~,i2] = max(corr);
end

