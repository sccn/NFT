% nft_dsl_inverse_problem_solution() - Distributed source localization
% (dsl) forward model generation
%
% Usage:
%   >> nft_dsl_inverse_problem_solution(subject_name, session_name, of)
%
% Inputs:
%   subject_name - subject name as entered in NFT main window
%   session_name - session name as entered in NFT main window
%   of - output folder where all the outputs are saved
%
% Optional keywords:
%
%   cond : conductivity vector (default = [0.33 0.0132 1.79 0.33])
%   mesh_name : mesh name that will be loaded (default: [subject_name '.1.msh'])
%   sensor_name:  sensor name (default: [subject_name '_' session_name '.sensors'])
%   ss_name : sourcespace name (default: [subject_name '_sourcespace.dip'])
%   LFM_name :  LFM name (default: session_name_LFM)
%
% Author: Zeynep Akalin Acar, SCCN, 2021

% Copyright (C) 2007 Zeynep Akalin Acar, SCCN, zeynep@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function nft_dsl_inverse_problem_solution(subject_name, session_name, of, EEG, comp_index, selection, sensor_file)

% selection = 2 -> SBL
% selection = 3 -> SCS 


% start source localization
curr_dir = pwd;
cd(of)


Phi_EEG = EEG.icawinv(:,comp_index);
ndip = length(comp_index);

% Solve distributed source loc for all models
sens = load(sensor_file, '-mat');
[ind_fp, ind_eeg] = find_indexes(EEG, sens.ind, sens.eloc);
Phi_EEG = Phi_EEG(ind_eeg,:);

LFM_name = [session_name '_LFM']; 
load(LFM_name)
LFM2 = LFM(ind_fp,:);
load ss_g10    
load FSss_cor    % sourcespace
load Node_area
max_scs_iter = 25;

if selection == 2
    % SBL
    load ss_g6    
    load ss_g3    
    disp('SBL source localization started...')
    for ij = 1:size(Phi_EEG, 2)
        Vdata = Phi_EEG(:,ij);
        Vdata = Vdata - mean(Vdata);

        Jb = source_loc_SBL_gaus_function(Vdata, ss_g3, ss_g6, ss_g10, LFM2, 4, 0.001, 0);
        sourceJ(:,ij) = Jb;
                
        pot = LFM2 * Jb; pot = pot - mean(pot);
        diff = Vdata - pot; 
        fvalJ(ij) = sum(diff(:).^2) / sum(Vdata(:).^2);

    end
    save cortex_source_sbl sourceJ fvalJ comp_index
    disp('SBL source localization finished...')
    
elseif selection == 3
    % SCS
    disp('SCS source localization started...')
    for ij = 1:size(Phi_EEG, 2)
        Vdata = Phi_EEG(:,ij);
        Vdata = Vdata - mean(Vdata);
        [Js1, Jit, fvx] = inverse_cov_sparse_average_noise_whole18_sparse_patchz2(LFM2, Vdata, ss_g10, max_scs_iter, 0);
        
        % find the most compact source in iterations
        for comi = 1:max_scs_iter+1
            pot = Jit(:,comi); pot = pot';
            compact0iter(comi) = calc_compactness(pot, An, Css(:,2:4), 1);
        end
        [maxcom, maxcomi] = max(compact0iter(5:max_scs_iter+1));
        maxcomiter = maxcomi + 4;

        sourceJ(:,ij) = Jit(:,maxcomiter);
                
        pot = LFM2 * Js1; pot = pot - mean(pot);
        diff = Vdata - pot; 
        fvalJ(ij) = sum(diff(:).^2) / sum(Vdata(:).^2);
    end
    save cortex_source_scs sourceJ fvalJ comp_index
    disp('SCS source localization finished...')
end

cd(curr_dir)


function [ind_fp, ind_eeg] = find_indexes(EEG, elp_index, eloc);
% realistic data icin

A = EEG.icawinv;

% neglect the FID electrodes
y = elp_index;
N = length(y);
for i = 1:N
    elocn(i).labels = eloc(y(i)).labels;
    elocn(i).X = eloc(y(i)).X;
    elocn(i).Y = eloc(y(i)).Y;
    elocn(i).Z = eloc(y(i)).Z;
 %   elocn(i).type = eloc(y(i)).type;
end

Nel2 = length(elocn);
Neeg = length(EEG.chanlocs);

Mel2 = zeros(1, Nel2);
Meeg = zeros(1, Neeg);

% Mel2(i) = index of EEG.chanlocs that correspond to electrode i of eloc2
% Meeg(i) = index of elocn that correspond to electrode i of EEG.chanlocs
clear Mel2 Meeg
for i = 1:Nel2
    for j = 1:Neeg
        if strcmp(elocn(i).labels, EEG.chanlocs(j).labels)
            Mel2(i) = j;
            Meeg(j) = i;
            continue;
        end
    end
end

ind_fp = find(Mel2>0); % index for the FP outputs (LFM, TM, session)
ind_eeg = find(Meeg>0); % index for the EEG structure (ICs)

function [J,Jit, fval,stdd_log_a,J_s,dispact_s,prob_ts] = inverse_cov_sparse_average_noise_whole18_sparse_patchz2(F,P,ss_MNI_gaussion,max_it,flag)

%J=inverse_cov_sparse_average(F,P,voxel_position)
%F is the lead fied matrix, P is the observed scalp potential and
% ss_MNI_gaussion for 6mm or 10mm
%max_it = 30, 20
% flag = 1
%voxel_position is the locatoin of the dipoles 
%Edited by Cheng Cao 2011
% A compact function is added 
% the covariance matrix is updated 
%Parallel computation is used 
% modified based on version 7, keep the hidden elements
%dealt with the noise issue, considering the DC shift of the noise
%Using two-point stepsize gradient
%Use log(std) to achieve better performance
%form version 10
%set the initial nsr_level according to the eigen value of M
%Under develovelpment
%Add smooth matrix Mar21,2012smooth_control
%change pinv to inv MAr 22

% initialize
[rt,ty] = size(F);
fval = zeros(1,max_it+1);
Jit = zeros(ty,max_it+1);

stop = 1;
step_size = 0.01; %0.01;
minium_nsr = 0.1; %0.1
n_control_para = 0.00001;
p_std_cof = 0;                              
nsr_coefi = 0.00005; %0.001;%0.01
%nsr_coefi = 0.1; %0.001;%0.01

[number_electrode,number_voxel] = size(F);
smooth_control = 0.1;
J_s = [];
dispact_s = [];
J = ones(number_voxel,1);
prob_ts = [];
MIN_GAMMA = 1e-16; 

P = P - mean(P);
P = P(1:number_electrode-1);
P_norm = norm(P);
P = P / P_norm;

F = F - ones(number_electrode,1) * mean(F);
F = F(1:number_electrode-1,:);
F_norms_sqr = (sum(F.^2))';

for i = 1:number_voxel
    F(:,i) = F(:,i)*((F_norms_sqr(i)).^-0.5);
end

F_e = zeros(number_electrode-1,number_electrode-1);
F_e = F_e-1 / number_electrode;


for i = 1:number_electrode-1
    F_e(i,i) = F_e(i,i)+1;
end

F_e = F_e / norm(F_e(:,1));

M = zeros(number_electrode-1,number_electrode-1);
N = zeros(number_voxel,1);
cov_column = zeros(number_voxel,1);
voxel_position_a = zeros(number_voxel,1);
pre_decompact = inf;

stdd_log_a = zeros(number_voxel,1);
nsr_level = zeros(number_electrode-1,1);
stdd_log_a = stdd_log_a+p_std_cof*log(F_norms_sqr)-0.1*log(min(F_norms_sqr));
n_it = 1;
err = zeros(1,number_electrode-1);

g_k = zeros(number_voxel+number_electrode-1,1);
g_k_old = g_k;
stdd_log_old = stdd_log_a;
nsr_level_old = nsr_level;
prob = 0;
   
J_index = 1:number_voxel;
number_a_voxel = number_voxel;
F_a = F;
n_itt = 1;

ss_matrix = ss_MNI_gaussion;
for i=1:number_voxel
    ss_MNI_gaussion(:,i) = ss_matrix(:,i) / sqrt(sum(ss_matrix(:,i).^2));
end
ss_MNI_gaussion = ss_MNI_gaussion';


while stop && n_it <= max_it
    
    row_temp = zeros(number_electrode-1,number_voxel);
    ss_diag_diag = sparse(1:number_voxel,1:number_voxel,exp(stdd_log_a));
    
    row_temp = F_a * ss_diag_diag * ss_MNI_gaussion;
    M = row_temp * row_temp';
    a1 = isnan(M);
    if sum(sum(a1))>1
        break
    end
    [~,D_m] = eig(M);
    if n_it == 1
        lam_max = max(diag(D_m));
        nsr_level_init = nsr_coefi*mean(diag(D_m));%Changed on Mar 19 2012
        scale = minium_nsr/nsr_level_init; 
      
        M = scale * M;
        row_temp = row_temp*sqrt(scale);
        stdd_log_a = stdd_log_a+0.5*log(scale);
        nsr_level = zeros(number_electrode-1,1)+log(minium_nsr)/2;
    end
    snr = mean(diag(D_m)) / exp(2*mean(nsr_level));
    M = M + F_e * diag(exp(2*nsr_level)) * F_e';
    [~,D_m]=eig(M);
    c_M = cond(M);
    min_eigv = min(diag(D_m));
    prob = P' * inv(M) * P;
    prob_t = (number_electrode-1) * log(P'*inv(M)*P) - log(det(inv(M))) - n_control_para * (mean(stdd_log_a) + mean(nsr_level));%%add the noise control;
    prob_ts = [prob_ts,prob_t];
    inv_M = inv(M);
   
    ss_diag_diag=sparse(1:number_voxel,1:number_voxel,exp(stdd_log_a));
   
    %%%%%%%%%%%%% calculate the current %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flag == 1
        N = F_a' * inv(M) * P;
        J_tmp = zeros(number_voxel,1);
        J_tmp = (ss_diag_diag*ss_MNI_gaussion)' * N;
        J = ss_diag_diag * ss_MNI_gaussion * J_tmp;
        rv = std(F*J-P)
        J = P_norm*J./sqrt(F_norms_sqr);

        if n_it > 1
            trange_J = norm(J-J_s(:,end)) / norm(J)
        end
        J_s = [J_s, J];
    end

    %%%%%%%%%%% starts to calcualte the gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    g_k = [];
    trace_inM = diag(inv_M);
    cov_init = (ss_MNI_gaussion);
    temp_PM=P'*inv_M;

    temp_PM_F=temp_PM*F;
    temp_cov_rowtemp_temp_PM=cov_init*row_temp'*temp_PM';
    g_1=-(number_electrode-1)*(2*temp_PM_F'.*temp_cov_rowtemp_temp_PM.*exp(stdd_log_a))/(prob);

    %temp_cov_rowtemp=cov_init*row_temp';
    temp_cov_rowtemp=row_temp*cov_init';
    g_2=sum(inv_M*F.*temp_cov_rowtemp);

    g_3=sum(inv_M*temp_cov_rowtemp.*F);
    g_k=[g_k;g_1+(g_3+g_2)'.*exp(stdd_log_a)];
    g_k=g_k+n_control_para/(number_voxel);

    for i=number_a_voxel+1:number_a_voxel+number_electrode-1
        ro_temp=zeros(1,number_electrode-1);
        ro_temp(i-number_a_voxel)=1;
        g_k(i)=-2*(number_electrode-1)*P'*inv_M*F_e*diag(ro_temp)*F_e'*inv_M*P*(exp(2*nsr_level(i-number_a_voxel)))/(prob)+2*trace(inv_M*F_e*diag(ro_temp)*F_e')*(exp(2*nsr_level(i-number_a_voxel)))-n_control_para/(number_electrode);
    end

    %fprintf('finished the calucation of the gradient\n');
    %%%%%%%%%%%%%%%%% transfer the gradient %%%%%%%%%%%%%%%%%%%%%%%%
    g_k_t=[];
    cov_init=ss_MNI_gaussion;
    g_k_t=[g_k_t;cov_init*g_k(1:number_voxel)];
    g_k_t=[g_k_t;g_k(number_a_voxel+1:end)];
    norm(g_k_t);
    surface_norm = ones(number_voxel,1);
    surface_norm = surface_norm/norm(surface_norm);

    g_k_t = g_k_t-mean(g_k_t);
    g_k = g_k'*g_k_t*g_k_t / (norm(g_k_t))^2;

    if n_itt >= 3
        step_size = min(0.5*sqrt(((stdd_log_a-stdd_log_old)'*(stdd_log_a-stdd_log_old)+(norm(nsr_level-nsr_level_old))^2)/((g_k-g_k_old)'*(g_k-g_k_old))),1000);
    end

    stdd_log_old = stdd_log_a;   
    nsr_level_old = nsr_level;

    stdd_log_a = stdd_log_a-step_size*g_k(1:number_voxel);
    nsr_level = nsr_level-step_size*g_k(number_voxel+1:end);

    M = zeros(number_electrode,number_electrode);

    g_k_old = g_k;
    n_it = n_it+1;
    n_itt = n_itt+1;
    %fprintf('%d iter  %d voxels used snr_level= %d current gradient',n_it-1,number_a_voxel,max(exp(nsr_level)),norm(g_k));
    std_max = max([stdd_log_a;nsr_level]);
    if std_max> 5 
        nsr_level = nsr_level-std_max+1;
        stdd_log_a = stdd_log_a-std_max+1;
        stdd_log_old = stdd_log_old-std_max+1;
        nsr_level_old = nsr_level_old-std_max+1;
        %fprintf('sum of std changed\n');
    end
    if norm(g_k) < 0.5   %||isinf(abs(prob_t))
        %stop = 0;
    end

   if norm(g_k) > 1000  % zeynep singular matrix oluyor
        stop = 0;
    end
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ss_diag_diag = sparse(1:number_voxel,1:number_voxel,exp(stdd_log_a));
    row_temp = zeros(number_electrode-1,number_voxel);

    cov_init = ss_MNI_gaussion;
    cov_init = ss_diag_diag*cov_init;

    row_temp = F_a * cov_init;
    M = row_temp * row_temp';
    M = M + F_e * diag(exp(2*nsr_level)) * F_e';
    N = F_a' * inv(M) * P;

    J_tmp = zeros(number_voxel,1);
    J_tmp = cov_init' * N;
    J = cov_init*J_tmp;

    err = F*J - P;
    fval(n_it) = sum(err(:).^2) / sum(P(:).^2);
    J = P_norm*J./sqrt(F_norms_sqr);
    Jit(:,n_it) = J;
    %if n_it > 5
    %    if fval(n_it)>fval(n_it-1)
    %        J = Jit(:,n_it-1);
    %        break;
    %    end
    %end

end
J_s = [J_s,J];


function compact = calc_compactness(pot, An, vr, flag, geodis)
% Cr, Er : cortical mesh
% An = area mesh nodes
% pot: cortical source
% vr: voxel_position

%flag = 0; % Cheng's method
%flag = 1; % Zunic's Kst
%flag = 2; % my area method
if nargin < 5
    geodis = [];
end

An  = An(:); An = An';
if flag == 0 
    k = 0.1;
    max_pot = max(abs(pot));
  
    ind_large_value = find(abs(pot) > max_pot/10);
     
    number_main_voxel = length(ind_large_value);
    if  number_main_voxel > 0.05*length(pot)
        compact = inf;
    else
        compact = 0;
        for i = 1:number_main_voxel-1
            for j = i+1:number_main_voxel
                te = exp(k*norm(vr(ind_large_value(j),:)-vr(ind_large_value(i),:)));
                compact = compact + abs(pot(ind_large_value(i))/max_pot) * abs(pot(ind_large_value(j))/max_pot) * te;
            end
        end
        compact = compact / (number_main_voxel*number_main_voxel);
    end
end



if flag == 3 
    k = 0.1;
    max_pot = max(abs(pot));
  
    ind_large_value = find(abs(pot) > max_pot/10);
     
    number_main_voxel = length(ind_large_value);
    if  number_main_voxel > 0.05*length(pot)
        compact = inf;
    else
        compact = 0;
        for i = 1:number_main_voxel-1
            for j = i+1:number_main_voxel
                if i > j
                    dist = geodis(i,j);
                else
                    dist = geodis(j,i);
                end
                if dist == 0
                    dist = 1000;
                end
                if i == j
                    dist = 0;
                end
                te = exp(k * dist);
                compact = compact + abs(pot(ind_large_value(i))/max_pot) * abs(pot(ind_large_value(j))/max_pot) * te;
            end
        end
        compact = compact / (number_main_voxel*number_main_voxel);
    end
end

if flag == 1
    % Compactness measure for 3D shapes
    % J Zunic, K. Hirota, C. Martinez-Ortiz, 2012
    pot = abs(pot);
      if max(pot) == 0
          compact = 0;
          return
      end
 
    pot = pot / max(pot) * 100; % max =100

    max_pot = max(abs(pot));
    bl = 10;
    ind_large_value = find(abs(pot) > max_pot/bl);
    while length(ind_large_value) < 2
        bl = bl*2;
        ind_large_value = find(abs(pot) > max_pot/bl);
    end
    number_main_voxel = length(ind_large_value);
    pot = abs(pot);
    %if  number_main_voxel > 0.05*length(pot)
        %compact = inf;
    %else
        vol_S = sum(An(ind_large_value).*pot(ind_large_value));
        surf_S = sum(An(ind_large_value));
        compact = 36 * pi * vol_S^2 / surf_S^3;
    %end
    compact = compact/max(pdist(vr(ind_large_value,:)));
end



if flag == 2
    pot = abs(pot);
    if max(pot) == 0
          compact = 0;
          return
      end
    pot = pot / max(pot) * 100; % max =100
    ap  = round(prctile(pot,99)); 
    if ap < 1; ap = 1; end
    ni = find(pot > ap);
    clear ai ao
    for i = ap:100
        
        ni = find(pot > i);
        ai(i) = sum(An(ni));
        ao(i) = sum(An(ni) .* pot(ni));
    end
    k = find(ao>std(ao)); jk = max(k)+1;
    compact = (1-ao(jk)/ao(ap)) * 100;
end
function [J_esmtime, Jt, fval]  = source_loc_SBL_gaus_function(p,ss_MNI_3_gaus,ss_MNI_6_gaus,ss_MNI_10_gaus,LFM_MNI_sourcespaceNn,max_inter,nsr_lv, flag1)

% if flag1 = 0, EM, if flag1 = 2; McKay

if nargin < 8
    flag1 = 0; % EM
end


LFM = LFM_MNI_sourcespaceNn;
%nsr_lv=0.0001;

ss_10 = sum(ss_MNI_10_gaus,2);
ss_6 = sum(ss_MNI_6_gaus,2);
ss_3 = sum(ss_MNI_3_gaus,2);
ss_10_spnorm = ss_MNI_10_gaus';
ss_6_spnorm = ss_MNI_6_gaus';
ss_3_spnorm = ss_MNI_3_gaus';
ii = length(ss);
for i = 1:ii;
    ss_10_spnorm(:,i) = ss_10_spnorm(:,i) / ss_10(i);
    ss_6_spnorm(:,i) = ss_6_spnorm(:,i) / ss_6(i);
    ss_3_spnorm(:,i) = ss_3_spnorm(:,i) / ss_3(i);
end
ss_10_spnorm = ss_10_spnorm';
ss_6_spnorm = ss_6_spnorm';
ss_3_spnorm = ss_3_spnorm';

New_lfm_10n = LFM * ss_10_spnorm;
New_lfm_6n = LFM * ss_6_spnorm;
New_lfm_3n = LFM * ss_3_spnorm;

%lfm2 = [New_lfm_10 New_lfm_6 New_lfm_3];
lfm3 = [New_lfm_10n New_lfm_6n New_lfm_3n];
ssx = [ss_10_spnorm' ss_6_spnorm' ss_3_spnorm'];

clear New_lfm_* ss_* LFM


[mut3,mu,dmu,kk,gamma,fval] = sparse_learning_ss(lfm3, p, nsr_lv, max_inter, flag1, 0, 0, 1); % EM

itern = size(mut3,1);
vn = mut3(itern,:);
J_esmtime = ssx * vn';
Jt = ssx * mut3'; 
   


function [mut, mu,dmu,k,gamma,fval] = sparse_learning_ss(Phi,T,lambda,iters,flag1,flag2,flag3, nofig)
% *************************************************************************
% 
% *** PURPOSE *** 
% Implements generalized versions of SBL and FOCUSS for learning sparse
% representations from possibly overcomplete dictionaries.
%
%
% *** USAGE ***
% [mu,dmu,k,gamma] = sparse_learning(Phi,T,lambda,iters,flag1,flag2,flag3);
%
%
% *** INPUTS ***
% Phi       = N X M dictionary
% T         = N X L data matrix
% lambda    = scalar trade-off parameter (balances sparsity and data fit)
% iters     = maximum number of iterations
%
% flag1     = 0: fast Mackay-based SBL update rules
% flag1     = 1: fast EM-based SBL update rule
% flag1     = 2: traditional (slow but sometimes better) EM-based SBL update rule
% flag1     = [3 p]: FOCUSS algorithm using the p-valued quasi-norm
%
% flag2     = 0: regular initialization (equivalent to min. norm solution)
% flag2     = gamma0: initialize with gamma = gamma0, (M X 1) vector
%
% flag3     = display flag; 1 = show output, 0 = supress output
%
% *** OUTPUTS ***
% mu        = M X L matrix of weight estimates
% dmu       = delta-mu at convergence
% k         = number of iterations used
% gamma     = M X 1 vector of hyperparameter values
%
%
% *************************************************************************
% Written by:  David Wipf, david.wipf@mrsc.ucsf.edu
% *************************************************************************
    

% *** Control parameters ***
MIN_GAMMA       = 1e-16;  % 1e-4
MIN_DMU         = 1e-12;  
MAX_ITERS       = iters;
DISPLAY_FLAG    = flag3;     % Set to zero for no runtime screen printouts


% *** Initializations ***
[N M] = size(Phi); 
[N L] = size(T);

if (~flag2)         gamma = ones(M,1);    
else                gamma = flag2;  end;   
 
keep_list = [1:M]';
m = length(keep_list);
mu = zeros(M,L);
dmu = -1;
k = 0;

iter=0; % zeynep
fig = 0;
if nargin < 8
   figure
else
    fig=1;
end

% *** Learning loop ***
while (1)
iter=iter+1; % zeynep

    % *** Prune things as hyperparameters go to zero ***
    if (min(gamma) < MIN_GAMMA )
		index = find(gamma > MIN_GAMMA);
		gamma = gamma(index);
		Phi = Phi(:,index);
		keep_list = keep_list(index);
        m = length(gamma);
     
        if (m == 0)   break;  end;
    end;
    
    
    % *** Compute new weights ***
    G = repmat(sqrt(gamma)',N,1);
    PhiG = Phi.*G; 
    [U,S,V] = svd(PhiG,'econ');
    
    [d1,d2] = size(S);
    if (d1 > 1)     diag_S = diag(S);  
    else            diag_S = S(1);      end;
    
    U_scaled = U(:,1:min(N,m)).*repmat((diag_S./(diag_S.^2 + lambda + 1e-16))',N,1);       
    Xi = G'.*(V*U_scaled'); 
        
    mu_old = mu;
    mu = Xi*T; 

    temp = zeros(M,L);
    if (m > 0) temp(keep_list,:) = mu;  end;
    mut(iter,:,:) = temp; % zeynep
    di = ceil(sqrt(iters));
    if fig==0
        subplot(di,di,iter); plot(mu); % zeynep
    end
    
    pot = Phi * mu;
    diff = T - pot;
    err = sum(diff(:).^2) / sum(T(:).^2);
    fval(iter) = err;
    
    % *** Update hyperparameters ***
    gamma_old = gamma;
    mu2_bar = sum(abs(mu).^2,2);
    
    if (flag1(1) == 0)
        % MacKay fixed-point SBL
        R_diag = real( (sum(Xi.'.*Phi)).' );
        te = L*R_diag;
        if min(te) < MIN_GAMMA  % zeynep
            ind_te = find(te<MIN_GAMMA);
            te(ind_te) = MIN_GAMMA;
        end
        gamma = mu2_bar./te;  
        
    elseif (flag1(1) == 1)
        % Fast EM SBL
        R_diag = real( (sum(Xi.'.*Phi)).' );
        gamma = sqrt( gamma.*real(mu2_bar./(L*R_diag)) ); 
        
    elseif (flag1(1) == 2)
        % Traditional EM SBL
        PhiGsqr = PhiG.*G;
        Sigma_w_diag = real( gamma - ( sum(Xi.'.*PhiGsqr) ).' );
        gamma = mu2_bar/L + Sigma_w_diag;
        
    else
        % FOCUSS
        p = flag1(2);
        gamma = (mu2_bar/L).^(1-p/2);
    end;
    
    
    
    % *** Check stopping conditions, etc. ***
  	k = k+1;   
    if (DISPLAY_FLAG) disp(['iters: ',num2str(k),'   num coeffs: ',num2str(m), ...
            '   gamma change: ',num2str(max(abs(gamma - gamma_old))), ...
            '   fval: ',num2str(err)]); end;    
    
    if (k >= MAX_ITERS) break;  end;
    
    % zeynep
    if iter>5
    if (abs(fval(iter-1)-fval(iter)) < 0.0001) break; end; % zeynep 6/11/15
    end
    %
    
	if (size(mu) == size(mu_old))
        dmu = max(max(abs(mu_old - mu)));
        if (dmu < MIN_DMU)  break;  end;
    end;
   
end;


% *** Expand weights, hyperparameters ***
temp = zeros(M,1);
if (m > 0) temp(keep_list,1) = gamma;  end;
gamma = temp;

temp = zeros(M,L);
if (m > 0) temp(keep_list,:) = mu;  end;
mu = temp;
   
return;


