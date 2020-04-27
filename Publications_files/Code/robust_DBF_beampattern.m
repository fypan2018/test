
%       自适应波束形成 最优约束矢量 SDP-LCMV
clc; clear all; close all;
j = sqrt(-1);
addpath('C:\Program Files\MATLAB\R2010a\personal\cvx');     % 添加 cvx 路径
N = 10;
K = 500;
lambda = 0.32;
d = lambda/2;
% -------------  噪声   ---------------------
noise = randn(N,K);
noise = sqrt(N/trace(noise*noise'/K))*noise; disp(['噪声功率为：',num2str(trace(noise*noise'/K/N))]);
% -------------  干扰   ---------------------
JNR = 30;
theta_j = 30/180*pi;                 % 夹角为入射方向与线阵法向的夹角
A_j = exp(j*2*pi*d/lambda*(0:N-1)'*sin(theta_j));
amp_j = randn(1,K)*10^(JNR/20);
jam = A_j*amp_j;
jam = jam*10^(JNR/20)/sqrt(trace(jam*jam'/K)/trace(noise*noise'/K));
disp(['干噪比 JNR = ',num2str(10*log10(trace(jam*jam'/K)/trace(noise*noise'/K)))]);

% --------------  信号  ---------------------
SNR = 10;
theta_s = 2/180*pi;
A_s = exp(j*2*pi*d/lambda*(0:N-1)'*sin(theta_s));
% amp_s = 10^(SNR/20);
% sig = A_s*amp_s; disp(['信噪比 SNR = ',num2str(10*log10(trace(sig*sig')/trace(noise*noise'/K)))]);
% data = jam + noise;
% data(:,round(K/2)) = data(:,round(K/2)) + sig;
%
amp_s = randn(1,K)*10^(SNR/20);
sig2 = A_s*amp_s;
sig2 = sig2*10^(SNR/20)/sqrt(trace(sig2*sig2'/K)/trace(noise*noise'/K));
data = sig2 + jam + noise;

Rx = data*data'/K;
% Rx = 10^(SNR/10)*A_s*A_s' + 10^(JNR/10)*A_j*A_j' + eye(N);      % 精确协方差矩阵

theta_e = 0/180*pi;
% theta_e = theta_s;

A_est = exp(j*2*pi*d/lambda*(0:N-1)'*sin(theta_e));

Rx = (Rx+Rx')/2;
Rx_dl = Rx + 5*eye(N);
%
w_mvdr_smi = Rx^(-1)*A_est/(A_est'*Rx^(-1)*A_est);
w_mvdr_smi_dl = Rx_dl^(-1)*A_est/(A_est'*Rx_dl^(-1)*A_est);
%
theta_c = theta_e + [-3,-1,0,1,3]/180*pi;
C = exp(j*2*pi*d/lambda*(0:N-1)'*sin(theta_c));
u = ones(length(theta_c),1);                                uu(1,:) = u;    % 固定的全1约束
w_lcmv_fix = Rx^(-1)*C*(C'*Rx^(-1)*C)^(-1)*u;
w_lcmv_fix_dl = Rx_dl^(-1)*C*(C'*Rx_dl^(-1)*C)^(-1)*u;
%
%              “幅相线性约束”
u = (A_est'*C)'; u = u/max(abs(u));                         uu(2,:) = u;
w_lcmv_mplc = Rx^(-1)*C*(C'*Rx^(-1)*C)^(-1)*u;
w_lcmv_mplc_dl = Rx_dl^(-1)*C*(C'*Rx_dl^(-1)*C)^(-1)*u;
%
%               “响应矢量优化——复数域实现”
C = exp(j*2*pi*d/lambda*(0:N-1)'*sin(theta_c));
% P = inv(C'*inv(Rx)*C);
% P = inv(C'*inv(Rx+100*eye(N))*C);
P = inv(C'*C);
P = (P + P')/2;
r_num = size(C,2);
for i = 1 : r_num
    E1(:,:,i) = zeros(r_num,r_num);
    E1(i,i,i) = 1;
end
cvx_begin quiet
    variable U(r_num,r_num) hermitian
    minimize(trace(P*U))
    subject to
        for i = 1 : r_num
            Ei = E1(:,:,i);
            trace(Ei*U) >= 1;
        end
        U == semidefinite(r_num);
cvx_end
[V1,D1] = eig(U);
D1 = abs(diag(D1));
u = V1(:,5);
w_lcmv_opt = Rx^(-1)*C*(C'*Rx^(-1)*C)^(-1)*u;
w_lcmv_opt_dl = Rx_dl^(-1)*C*(C'*Rx_dl^(-1)*C)^(-1)*u;
%
%           “响应矢量优化——实数域实现”
P_ext = [real(P),-imag(P);imag(P),real(P)];
B = [eye(r_num),eye(r_num)];
cvx_begin quiet
    variable U(2*r_num,2*r_num) symmetric
    minimize(trace(P_ext*U))
    subject to
        for i = 1 : r_num
            Ei = E1(:,:,i);
            trace(B.'*Ei*B*U) >= 1;
        end
        U == semidefinite(2*r_num);
cvx_end
[V1,D1] = eig(U);
D1 = abs(diag(D1));
u = V1(:,2*r_num);
u = u(1:r_num) + j*u(r_num+1:end);
u = u/max(abs(u));                  uu(3,:) = u;    uu(4,:) = -u;   % 半正定规划最优约束
w_lcmv_opt2 = Rx^(-1)*C*(C'*Rx^(-1)*C)^(-1)*u;
w_lcmv_opt2_dl = Rx_dl^(-1)*C*(C'*Rx_dl^(-1)*C)^(-1)*u;

%          “响应矢量优化——低秩逼近实现”
% 当存在目标信号时，不稳健，这个怎么解释？如何利用 Samuel 的文章思想？
theta_c = theta_e + linspace(-6,6,9)/180*pi;
C = exp(j*2*pi*d/lambda*(0:N-1)'*sin(theta_c));
[Y,S,Z] = svd(C);
num_dof = 3;
Y = Y(:,1:num_dof);
S = S(1:num_dof,1:num_dof);
Z = Z(:,1:num_dof);

% P = Z*S^(-1)*inv(Y'*inv(Rx)*Y)*S^(-1)*Z';
P = Z*S^(-1)*inv(Y'*Y)*S^(-1)*Z';
P_ext = [real(P),-imag(P);imag(P),real(P)];
r_num = size(C,2);
for i = 1 : r_num
    E3(:,:,i) = zeros(r_num,r_num);
    E3(i,i,i) = 1;
end
B = [eye(r_num),eye(r_num)];
cvx_begin quiet
    variable U(2*r_num,2*r_num) symmetric
    minimize(trace(P_ext*U))
    subject to
        for i = 1 : r_num
            Ei = E3(:,:,i);
            trace(B.'*Ei*B*U) >= 1;
        end
        U == semidefinite(2*r_num);
cvx_end
[V1,D1] = eig(U);
D1 = abs(diag(D1));
u = V1(:,2*r_num);      % figure(); plot(10*log10(abs(D1)),'-mo'); hold on;
u = u(1:r_num) + j*u(r_num+1:end);
w_lcmv_opt3 = Rx^(-1)*Y*(Y'*Rx^(-1)*Y)^(-1)*S^(-1)*Z'*u;
%
%        “worst case”     by Zhiquan Luo
delta_fs = 1/2/N;
S_delta = exp(j*2*pi*(0:N-1)'*delta_fs);
epsilon_s = sqrt(norm(S_delta))*0.9;
% 优化函数
cvx_begin quiet
    variable u(N) complex
    minimize(u'*Rx*u)
    subject to
        epsilon_s*norm(u) <= real(A_est'*u)-1;
        imag(A_est'*u) == 0;
cvx_end
w_worst_case = u;
%
%        “RAB-SDP”    by Zhuliang Yu
theta_yu = theta_e + linspace(-3,3,100)/180*pi;
rdb = 0.3;
Upper = 10^(rdb/20);
Lower = 10^(-rdb/20);
cvx_begin quiet
    variable W(N,N) hermitian
    minimize(trace(Rx*W))
    subject to
        for i = 1 : length(theta_yu)
            A_s = exp(j*2*pi*d/lambda*(0:N-1)'*sin(theta_yu(i)));
            trace(A_s*A_s'*W) >= Lower^2;
            trace(A_s*A_s'*W) <= Upper^2;
        end
        W == semidefinite(N);
cvx_end

r = zeros(1,2*N-1);
for m = 1 : N
    for n = 1 : N
        index = m - n + N;
        r(index) = r(index) + W(m,n);
    end
end
num = 10*N;
T = zeros(num,num);
for m = 1 : num
    for n = 1 : num
        index = m-n + N;
        if index > 0 & index < 2*N
            T(m,n) = r(index);
        else
            T(m,n) = 0;
        end
    end
end
T=(T+T')/2;
R_comp = chol(T);    L_comp = R_comp';
w_SDP_yu1 = fliplr(L_comp(end,:));
w_SDP_yu1(N+1:num) = [];    w_SDP_yu1=w_SDP_yu1.'/norm(w_SDP_yu1);

[V1,D1] = eig(W);
w_SDP_yu2 = V1(:,N);
% [v,h] = eig(W,'nobalance');   h1 = sort(diag(h)); 
% figure(); xlabel('number of eigenvalue'); hold on;
% plot(10*log10((abs(h1))),'b+-');axis tight;

%
%           
theta = asin(linspace(-1,1,361));
for n = 1:length(theta)
    a = exp(j*2*pi*d/lambda*(0:N-1)'*sin(theta(n)));
    a = a/norm(a);
    Q1(n) = w_mvdr_smi'*(a*a')*w_mvdr_smi/norm(w_mvdr_smi)^2;
    Q2(n) = w_mvdr_smi_dl'*(a*a')*w_mvdr_smi_dl/norm(w_mvdr_smi_dl)^2;
    Q3(n) = w_lcmv_fix'*(a*a')*w_lcmv_fix/norm(w_lcmv_fix)^2;
    Q4(n) = w_lcmv_fix_dl'*(a*a')*w_lcmv_fix_dl/norm(w_lcmv_fix_dl)^2;
    Q5(n) = w_lcmv_mplc'*(a*a')*w_lcmv_mplc/norm(w_lcmv_mplc)^2;
    Q6(n) = w_lcmv_mplc_dl'*(a*a')*w_lcmv_mplc_dl/norm(w_lcmv_mplc_dl)^2;
    Q7(n) = w_lcmv_opt'*(a*a')*w_lcmv_opt/norm(w_lcmv_opt)^2;
    Q8(n) = w_lcmv_opt_dl'*(a*a')*w_lcmv_opt_dl/norm(w_lcmv_opt_dl)^2;
    Q9(n) = w_lcmv_opt2'*(a*a')*w_lcmv_opt2/norm(w_lcmv_opt2)^2;
    Q10(n) = w_lcmv_opt2_dl'*(a*a')*w_lcmv_opt2_dl/norm(w_lcmv_opt2_dl)^2;
    Q11(n) = w_lcmv_opt3'*(a*a')*w_lcmv_opt3/norm(w_lcmv_opt3)^2;
    Q12(n) = w_worst_case'*(a*a')*w_worst_case/norm(w_worst_case)^2;
    Q13(n) = w_SDP_yu1'*(a*a')*w_SDP_yu1/norm(w_SDP_yu1)^2;  
    Q14(n) = w_SDP_yu2'*(a*a')*w_SDP_yu2/norm(w_SDP_yu2)^2; 
end

Q1 = Q1/max(abs(Q1));
Q2 = Q2/max(abs(Q2));
Q3 = Q3/max(abs(Q3));
Q4 = Q4/max(abs(Q4));
Q5 = Q5/max(abs(Q5));
Q6 = Q6/max(abs(Q6));
Q7 = Q7/max(abs(Q7));
Q8 = Q8/max(abs(Q8));
Q9 = Q9/max(abs(Q9));
Q10 = Q10/max(abs(Q10));
Q11 = Q11/max(abs(Q11));
Q12 = Q12/max(abs(Q12));
Q13 = Q13/max(abs(Q13));
Q14 = Q14/max(abs(Q14));

m = find(min(abs(theta-theta_e)) == abs(theta-theta_e));
Q1 = Q1/abs(Q1(m));
Q2 = Q2/max(abs(Q2));
Q3 = Q3/abs(Q3(m));
Q4 = Q4/abs(Q4(m));
Q5 = Q5/abs(Q5(m));
Q6 = Q6/abs(Q6(m));
Q7 = Q7/abs(Q7(m));
Q8 = Q8/abs(Q8(m));
Q9 = Q9/abs(Q9(m));
Q10 = Q10/abs(Q10(m));
Q11 = Q11/abs(Q11(m));
Q12 = Q12/abs(Q12(m));
Q13 = Q13/abs(Q13(m));
Q14 = Q14/abs(Q14(m));

figure(); xlabel('DOA(deg)');ylabel('Beampattern (dB)'); hold on; box on;
plot(theta*180/pi,10*log10(abs(Q1)),':g','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q2)),'-g','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q3)),':b','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q4)),'-b','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q5)),':k','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q6)),'-k','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q7)),'--m','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q8)),'-m','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q9)),'--r','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q10)),'-r','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q11)),'-.r','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q12)),'-.c','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q13)),'--c','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q14)),'-c','LineWidth',2);

legend('MVDR-SMI','MVDR-LSMI','LCMV','DL-LCMV','MPLC LCMV','MPLC LCMV','LCMV-SDP','LCMV-SDP','LCMV-SDP','LCMV-SDP','LCMV-SDP low rank','worst case','RAB-SDP','RAB-SDP',3);
plot(theta_s/pi*180*ones(1,100),linspace(min(10*log10(abs(Q1))),max(10*log10(abs(Q1)))+10,100),'--m');hold on;    %信号真实来向
plot(theta_j/pi*180*ones(1,100),linspace(min(10*log10(abs(Q1))),max(10*log10(abs(Q1)))+10,100),'--m');hold on; %干扰1方向


%%             文章用图   方向图仿真
figure(); xlabel('DOA(deg)');ylabel('Beampattern (dB)'); hold on; box on;
plot(theta*180/pi,10*log10(abs(Q1)),':g','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q12)),'-.c','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q3)),'--k','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q4)),'-b','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q9)),'--.m','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q10)),'-..r','LineWidth',2);
legend('MVDR-SMI','worst case','LCMV','DL LCMV','SDR LCMV','DL SDR LCMV',3);
plot(theta_s/pi*180*ones(1,100),linspace(min(10*log10(abs(Q1))),max(10*log10(abs(Q1)))+10,100),':k');hold on;    %信号真实来向
plot(theta_j/pi*180*ones(1,100),linspace(min(10*log10(abs(Q1))),max(10*log10(abs(Q1)))+10,100),':k');hold on; %干扰1方向
axis([-100,100,-80,10]);

figure(); xlabel('DOA(deg)');ylabel('Beampattern (dB)'); hold on; box on;
plot(theta*180/pi,10*log10(abs(Q2)),':r','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q4)),'--b','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q12)),'-.g','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q13)),'-k','LineWidth',2);
plot(theta*180/pi,10*log10(abs(Q10)),'--.m','LineWidth',2);
legend('LSMI-MVDR','Traditional LCMV','Worst case optimization','RAB-SDP','Proposed method',3);
plot(theta_s/pi*180*ones(1,100),linspace(min(10*log10(abs(Q1))),max(10*log10(abs(Q1)))+10,100),':k');hold on;    %信号真实来向
plot(theta_j/pi*180*ones(1,100),linspace(min(10*log10(abs(Q1))),max(10*log10(abs(Q1)))+10,100),':k');hold on; %干扰1方向
axis([-100,100,-80,10]);
