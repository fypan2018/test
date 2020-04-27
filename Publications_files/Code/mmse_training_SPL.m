clear all,
close all,

K    = 6;                   % Number of digital antenna channels
snr  = 0;
theta    = d2r([0 5.5]);    %  Target angles. The number of targets is
Sample   =  [20 30 40 60 80 100 120 150 200 300 400 600 900 1200]';  % Number of snapshots
C = 2;
Ntrial  = 100;

jj = sqrt(-1);
% ----------------------------------------------------------------------- %
% Commands.
% ----------------------------------------------------------------------- %
Ndoa         = size(theta,2);
Nsample      = length(Sample);
pdf_MDL        = zeros(Nsample,Ndoa+2);
proposedMethod = zeros(Nsample,Ndoa+2);
Num_ref      = [0:Ndoa+1]';
for nNsample =1:Nsample
    number_dEVD     = zeros(Ndoa+2,1);
    number_CS       = zeros(Ndoa+2,1);

    SNR = ones(size(theta',1),1) * snr;
    T = Sample(nNsample);
    for nTrial = 1:Ntrial
        %==================================================================
        source_power = (10.^(SNR./10));
        source_amplitude = sqrt(source_power)*ones(1,T); % Source standard deviation.
        source_wave = sqrt(0.5)*(randn(T,Ndoa) + jj*randn(T,Ndoa));
        st = source_amplitude.*source_wave.';
        d0 = st(1,:).';
        nt = sqrt(0.5)*(randn(K,T)+jj*randn(K,T));
        A = exp(jj*pi*[0:K-1]'*sin(theta));
        xt = A*st + nt;
        %============================== MMDL method ========================
        xi = C*sqrt((log(log(T)))/T);
        Rxx = xt*xt'./T;
        rxd = xt*conj(d0)./T;
        [ev,dd,eu] = svd(Rxx);
        dd = (diag(dd));
        cv = abs(ev'*rxd);
        rMM = abs(d0'*d0)./T;
        dk = (cv.^2./dd) ;
        Detector_MDL_Gdisk = zeros(K,1);
        %------------------------------------------
        for k = 1:K
            beta(k) = sum(dk(1:k));
            MMSE0(k) = (rMM - beta(k));
        end
        Test0 = abs(mean(MMSE0) - rMM);
        if Test0 < xi
            Enumber_CS = 0;
        else
            for k = 1:K
                MMSE(k) = T*log(MMSE0(k));
                Pay(k) = (k^2+k)*log(T)./2;
                Detector_MDL_Gdisk(k) = MMSE(k) + Pay(k);
            end
            [Evalue, Enumber_CS] = min(Detector_MDL_Gdisk);
        end
        %------------------------------------------
        p_CS = find(Enumber_CS - Num_ref == 0);
        number_CS(p_CS) = number_CS(p_CS) + 1;
        %======================= CMDL method ==============================
        [Ke,N]=size(xt);
        Rx = (xt*xt')./T;
        [u,s,v] = svd(Rx);
        sd = diag(s);
        a = zeros(1,K);
        for m = 0:K-1
            negv = sd(m+1:K);
            Tsph = mean(negv)/((prod(negv))^(1/(Ke-m)));
            a(m+1) = T*(K-m)*log(Tsph) + m*(2*K-m)*log(T)/2;
        end
        [y,b] = min(a);
        dEVD = b - 1;

        p_dEVD = find(dEVD - Num_ref == 0);
        number_dEVD(p_dEVD) = number_dEVD(p_dEVD) + 1;
        %----------------------------------------------
    end %for nTrial
    pdf_MDL(nNsample,1:end)         = number_dEVD'/nTrial;
    proposedMethod(nNsample,1:end)  = number_CS'/nTrial;
end %
%============================================
figure,
semilogx(Sample,pdf_MDL(:,Ndoa+1),'b:*',Sample,proposedMethod(:,Ndoa+1),'r-o')
legend('CMDL','Proposed MMDL',4)
ylabel('Probability of Detection')
xlabel('Number of Snapshots')
axis([Sample(1),Sample(length(Sample)),0,1])
