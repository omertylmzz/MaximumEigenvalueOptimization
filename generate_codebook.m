function codebook = generate_codebook(sim_parameters)
    phi = @(x)exp(1i*pi*x/2);

    XP = sim_parameters.nof_ant_ports_XP;
    N1 = sim_parameters.nof_ant_ports_N1;
    N2 = sim_parameters.nof_ant_ports_N2;
    codebook_mode = sim_parameters.codebook_mode;
    nof_ant_ports = XP * N1 * N2;
    nof_layers = sim_parameters.nof_layers;
    O1 = 4;
    if N2 > 1; O2 = 4; else; O2 = 1; end

    if nof_ant_ports == 2
        if nof_layers == 1
            codebook(:,:,1) = 1/sqrt(2).*[1; 1];
            codebook(:,:,2) = 1/sqrt(2).*[1; 1i];
            codebook(:,:,3) = 1/sqrt(2).*[1; -1];
            codebook(:,:,4) = 1/sqrt(2).*[1; -1i];
        elseif nof_layers == 2
            codebook(:,:,1) = 1/2*[1 1;1 -1];
            codebook(:,:,2) = 1/2*[1 1; 1i -1i];
        end
    elseif nof_ant_ports > 2
        switch nof_layers
            case 1 
                i11_length = N1*O1;
                i12_length = N2*O2;
                i2_length = 4;
                codebook = zeros(nof_ant_ports,nof_layers,i2_length,i11_length,i12_length);
                for i11 = 0:i11_length-1
                    for i12 = 0:i12_length-1
                        for i2 = 0:i2_length-1
                            l = i11;
                            m = i12;
                            n = i2;
                            vlm = get_vlm(N1,N2,O1,O2,l,m);
                            phi_n = phi(n);
                            codebook(:,:,i2+1,i11+1,i12+1) = (1/sqrt(nof_ant_ports))*[vlm ;...
                                phi_n*vlm];
                        end
                    end
                end

            case 2 
                if (N1 > N2) && (N2 > 1)
                    i13_length = 4;
                    k1 = [0 O1 0 2*O1];
                    k2 = [0 0 O2 0];
                elseif N1 == N2
                    i13_length = 4;
                    k1 = [0 O1 0 O1];
                    k2 = [0 0 O2 O2];
                elseif (N1 == 2) && (N2 == 1)
                    i13_length = 2;
                    k1 = O1*(0:1);
                    k2 = [0 0];
                else
                    i13_length = 4;
                    k1 = O1*(0:3);
                    k2 = [0 0 0 0] ;
                end

                if codebook_mode == 1
                    i11_length = N1*O1;
                    i12_length = N2*O2;
                    i2_length = 2;
                    codebook = zeros(nof_ant_ports,nof_layers,i2_length,i11_length,i12_length,i13_length);
                    for i11 = 0:i11_length-1
                        for i12 = 0:i12_length-1
                            for i13 = 0:i13_length-1
                                for i2 = 0:i2_length-1
                                    l = i11;
                                    m = i12;
                                    n = i2;
                                    lPrime = i11+k1(i13+1);
                                    mPrime = i12+k2(i13+1);
                                    vlm = get_vlm(N1,N2,O1,O2,l,m);
                                    vlPrime_mPrime = get_vlm(N1,N2,O1,O2,lPrime,mPrime);
                                    phi_n = phi(n);
                                    codebook(:,:,i2+1,i11+1,i12+1,i13+1) = ...
                                        (1/sqrt(2*nof_ant_ports))*[vlm        vlPrime_mPrime;...
                                        phi_n*vlm  -phi_n*vlPrime_mPrime];
                                end
                            end
                        end
                    end
                else 
                    i11_length = N1*O1/2;
                    if N2 == 1
                        i12_length = 1;
                    else
                        i12_length = N2*O2/2;
                    end
                    i2_length = 8;
                    codebook = zeros(nof_ant_ports,nof_layers,i2_length,i11_length,i12_length,i13_length);
                    for i11 = 0:i11_length-1
                        for i12 = 0:i12_length-1
                            for i13 = 0:i13_length-1
                                for i2 = 0:i2_length-1
                                    floor_i2by2 = floor(i2/2);
                                    if N2 == 1
                                        l = 2*i11 + floor_i2by2;
                                        lPrime = 2*i11 + floor_i2by2 + k1(i13+1);
                                        m = 0;
                                        mPrime = 0;
                                    else % N2 > 1
                                        lmAddVals = [0 0; 1 0; 0 1;1 1];
                                        l = 2*i11 + lmAddVals(floor_i2by2+1,1);
                                        lPrime =  2*i11 + k1(i13+1) + lmAddVals(floor_i2by2+1,1);
                                        m = 2*i12 + lmAddVals(floor_i2by2+1,2);
                                        mPrime =  2*i12 + k2(i13+1) + lmAddVals(floor_i2by2+1,2);
                                    end
                                    n = mod(i2,2);
                                    vlm = get_vlm(N1,N2,O1,O2,l,m);
                                    vlPrime_mPrime = get_vlm(N1,N2,O1,O2,lPrime,mPrime);
                                    phi_n = phi(n);
                                    codebook(:,:,i2+1,i11+1,i12+1,i13+1) = ...
                                        (1/sqrt(2*nof_ant_ports))*[vlm        vlPrime_mPrime;...
                                        phi_n*vlm  -phi_n*vlPrime_mPrime];
                                end
                            end
                        end
                    end
                end

            case {3,4} 
                if (nof_ant_ports < 16)
                    if (N1 == 2) && (N2 == 1)
                        i13_length = 1;
                        k1 = O1;
                        k2 = 0;
                    elseif (N1 == 4) && (N2 == 1)
                        i13_length = 3;
                        k1 = O1*(1:3);
                        k2 = [0 0 0];
                    elseif (N1 == 6) && (N2 == 1)
                        i13_length = 4;
                        k1 = O1*(1:4);
                        k2 = [0 0 0 0];
                    elseif (N1 == 2) && (N2 == 2)
                        i13_length = 3;
                        k1 = [O1 0 O1];
                        k2 = [0 O2 O2];
                    elseif (N1 == 3) && (N2 == 2)
                        i13_length = 4;
                        k1 = [O1 0 O1 2*O1];
                        k2 = [0 O2 O2 0];
                    end

                    i11_length = N1*O1;
                    i12_length = N2*O2;
                    i2_length = 2;
                    codebook = zeros(nof_ant_ports,nof_layers,i2_length,i11_length,i12_length,i13_length);
                    for i11 = 0:i11_length-1
                        for i12 = 0:i12_length-1
                            for i13 = 0:i13_length-1
                                for i2 = 0:i2_length-1
                                    l = i11;
                                    lPrime = i11+k1(i13+1);
                                    m = i12;
                                    mPrime = i12+k2(i13+1);
                                    n = i2;
                                    vlm = get_vlm(N1,N2,O1,O2,l,m);
                                    vlPrime_mPrime = get_vlm(N1,N2,O1,O2,lPrime,mPrime);
                                    phi_n = phi(n);
                                    phi_vlm = phi_n*vlm;
                                    phi_vlPrime_mPrime = phi_n*vlPrime_mPrime;
                                    if nof_layers == 3
                                        codebook(:,:,i2+1,i11+1,i12+1,i13+1) = ...
                                            (1/sqrt(3*nof_ant_ports))*[vlm      vlPrime_mPrime      vlm;...
                                            phi_vlm  phi_vlPrime_mPrime  -phi_vlm];
                                    else
                                        codebook(:,:,i2+1,i11+1,i12+1,i13+1) = ...
                                            (1/sqrt(4*nof_ant_ports))*[vlm      vlPrime_mPrime      vlm       vlPrime_mPrime;...
                                            phi_vlm  phi_vlPrime_mPrime  -phi_vlm  -phi_vlPrime_mPrime];
                                    end
                                end
                            end
                        end
                    end
                else 
                    i11_length = N1*O1/2;
                    i12_length = N2*O2;
                    i13_length = 4;
                    i2_length = 2;
                    codebook = zeros(nof_ant_ports,nof_layers,i2_length,i11_length,i12_length,i13_length);
                    for i11 = 0:i11_length-1
                        for i12 = 0:i12_length-1
                            for i13 = 0:i13_length-1
                                for i2 = 0:i2_length-1
                                    theta = exp(1i*pi*i13/4);
                                    l = i11;
                                    m = i12;
                                    n = i2;
                                    phi_n = phi(n);
                                    vbarlm = get_vbarlm(N1,N2,O1,O2,l,m);
                                    theta_vbarlm = theta*vbarlm;
                                    phi_vbarlm = phi_n*vbarlm;
                                    phi_theta_vbarlm = phi_n*theta*vbarlm;
                                    if nof_layers == 3
                                        codebook(:,:,i2+1,i11+1,i12+1,i13+1) = ...
                                            (1/sqrt(3*nof_ant_ports))*[vbarlm            vbarlm             vbarlm;...
                                            theta_vbarlm      -theta_vbarlm      theta_vbarlm;...
                                            phi_vbarlm        phi_vbarlm         -phi_vbarlm;...
                                            phi_theta_vbarlm  -phi_theta_vbarlm  -phi_theta_vbarlm];
                                    else
                                        codebook(:,:,i2+1,i11+1,i12+1,i13+1) = ...
                                            (1/sqrt(4*nof_ant_ports))*[vbarlm            vbarlm             vbarlm             vbarlm;...
                                            theta_vbarlm      -theta_vbarlm      theta_vbarlm       -theta_vbarlm;...
                                            phi_vbarlm        phi_vbarlm         -phi_vbarlm        -phi_vbarlm;...
                                            phi_theta_vbarlm  -phi_theta_vbarlm  -phi_theta_vbarlm  phi_theta_vbarlm];
                                    end
                                end
                            end
                        end
                    end
                end

            case {5,6} 
                i11_length = N1*O1;
                if N2 == 1
                    i12_length = 1;
                else % N2 > 1
                    i12_length = N2*O2;
                end
                i2_length = 2;
                codebook = zeros(nof_ant_ports,nof_layers,i2_length,i11_length,i12_length);
                for i11 = 0:i11_length-1
                    for i12 = 0:i12_length-1
                        for i2 = 0:i2_length-1
                            if N2 == 1
                                l = i11;
                                lPrime = i11+O1;
                                l_dPrime = i11+2*O1;
                                m = 0;
                                mPrime = 0;
                                m_dPrime = 0;
                            else % N2 > 1
                                l = i11;
                                lPrime = i11+O1;
                                l_dPrime = i11+O1;
                                m = i12;
                                mPrime = i12;
                                m_dPrime = i12+O2;
                            end
                            n = i2;
                            vlm = get_vlm(N1,N2,O1,O2,l,m);
                            vlPrime_mPrime = get_vlm(N1,N2,O1,O2,lPrime,mPrime);
                            vlDPrime_mDPrime = get_vlm(N1,N2,O1,O2,l_dPrime,m_dPrime);
                            phi_n = phi(n);
                            phi_vlm = phi_n*vlm;
                            phi_vlPrime_mPrime = phi_n*vlPrime_mPrime;
                            if nof_layers == 5
                                codebook(:,:,i2+1,i11+1,i12+1) = ...
                                    1/(sqrt(5*nof_ant_ports))*[vlm       vlm        vlPrime_mPrime   vlPrime_mPrime    vlDPrime_mDPrime;...
                                    phi_vlm   -phi_vlm   vlPrime_mPrime   -vlPrime_mPrime   vlDPrime_mDPrime];
                            else
                                codebook(:,:,i2+1,i11+1,i12+1) = ...
                                    1/(sqrt(6*nof_ant_ports))*[vlm       vlm        vlPrime_mPrime       vlPrime_mPrime        vlDPrime_mDPrime   vlDPrime_mDPrime;...
                                    phi_vlm   -phi_vlm   phi_vlPrime_mPrime   -phi_vlPrime_mPrime   vlDPrime_mDPrime   -vlDPrime_mDPrime];
                            end
                        end
                    end
                end

            case{7,8} 
                if N2 == 1
                    i12_length = 1;
                    if N1 == 4
                        i11_length = N1*O1/2;
                    else % N1 > 4
                        i11_length = N1*O1;
                    end
                else % N2 > 1
                    i11_length = N1*O1;
                    if (N1 == 2 && N2 == 2) || (N1 > 2 && N2 > 2)
                        i12_length = N2*O2;
                    else % (N1 > 2 && N2 == 2)
                        i12_length = N2*O2/2;
                    end
                end
                i2_length = 2;
                codebook = zeros(nof_ant_ports,nof_layers,i2_length,i11_length,i12_length);
                for i11 = 0:i11_length-1
                    for i12 = 0:i12_length-1
                        for i2 = 0:i2_length-1
                            if N2 == 1
                                l = i11;
                                lPrime = i11+O1;
                                l_dPrime = i11+2*O1;
                                l_tPrime = i11+3*O1;
                                m = 0;
                                mPrime = 0;
                                m_dPrime = 0;
                                m_tPrime = 0;
                            else % N2 > 1
                                l = i11;
                                lPrime = i11+O1;
                                l_dPrime = i11;
                                l_tPrime = i11+O1;
                                m = i12;
                                mPrime = i12;
                                m_dPrime = i12+O2;
                                m_tPrime = i12+O2;
                            end
                            n = i2;
                            vlm = get_vlm(N1,N2,O1,O2,l,m);
                            vlPrime_mPrime = get_vlm(N1,N2,O1,O2,lPrime,mPrime);
                            vlDPrime_mDPrime = get_vlm(N1,N2,O1,O2,l_dPrime,m_dPrime);
                            vlTPrime_mTPrime = get_vlm(N1,N2,O1,O2,l_tPrime,m_tPrime);
                            phi_n = phi(n);
                            phi_vlm = phi_n*vlm;
                            phi_vlPrime_mPrime = phi_n*vlPrime_mPrime;
                            if nof_layers == 7
                                codebook(:,:,i2+1,i11+1,i12+1) = ...
                                    1/(sqrt(7*nof_ant_ports))*[vlm       vlm        vlPrime_mPrime       vlDPrime_mDPrime   vlDPrime_mDPrime    vlTPrime_mTPrime   vlTPrime_mTPrime;...
                                    phi_vlm   -phi_vlm   phi_vlPrime_mPrime   vlDPrime_mDPrime   -vlDPrime_mDPrime   vlTPrime_mTPrime   -vlTPrime_mTPrime];
                            else
                                codebook(:,:,i2+1,i11+1,i12+1) = ...
                                    1/(sqrt(8*nof_ant_ports))*[vlm       vlm        vlPrime_mPrime       vlPrime_mPrime        vlDPrime_mDPrime   vlDPrime_mDPrime    vlTPrime_mTPrime   vlTPrime_mTPrime;...
                                    phi_vlm   -phi_vlm   phi_vlPrime_mPrime   -phi_vlPrime_mPrime   vlDPrime_mDPrime   -vlDPrime_mDPrime   vlTPrime_mTPrime   -vlTPrime_mTPrime];
                            end
                        end
                    end
                end
        end
    end    
end

function vlm = get_vlm(N1,N2,O1,O2,l,m)
    um = exp(2*pi*1i*m*(0:N2-1)/(O2*N2));
    ul = exp(2*pi*1i*l*(0:N1-1)/(O1*N1)).';
    vlm =  reshape((ul.*um).',[],1);
end

function vbarlm = get_vbarlm(N1,N2,O1,O2,l,m)
    um = exp(2*pi*1i*m*(0:N2-1)/(O2*N2));
    ul = exp(2*pi*1i*l*(0:(N1/2)-1)/(O1*N1/2)).';
    vbarlm = reshape((ul.*um).',[],1);
end
