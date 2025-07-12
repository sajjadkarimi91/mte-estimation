
function [PTE ,I , PMI , NCS , sig_te]...
    = te_lsim( pi_0_lsim , coupling_tetha_IM , transition_matrices_IM ,  gmm_para_lsim  ,  Y_in, extra  )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Transfer Entropy estimation via LSIMs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C = size(Y_in,1);
flag_exact = extra.flag_exact;
min_order = extra.min_order;

Y_in_temp = Y_in;
for c = 1:C
    Y_in_temp{c,1} = nan*Y_in_temp{c,1};
end

[~ ,alpha_out , ~ , ~ , ~]...
    = mc_lsim( pi_0_lsim , coupling_tetha_IM , transition_matrices_IM ,  gmm_para_lsim  ,  Y_in_temp  );

for c = 1:C
    pi_steady_lsim{c,1} = alpha_out{c,1}(:,end);
end

[P_ot_c_cond_past_temp ,alpha_steady ]...
    = mc_lsim( pi_steady_lsim , coupling_tetha_IM , transition_matrices_IM ,  gmm_para_lsim  ,  Y_in_temp  );

[P_ot_c_cond_past ,alpha , ~ , P_observ_cond_to_state_out ]...
    = mc_lsim( pi_steady_lsim , coupling_tetha_IM , transition_matrices_IM ,  gmm_para_lsim  ,  Y_in  );

H_ot_Opast = - mean(log(P_ot_c_cond_past(:,min_order:end)),2,'omitnan');

for c = 1:C
    P_ot_c_cond_past_temp(c,:) = sum(alpha_steady{c,1}.*P_observ_cond_to_state_out{c,1});
    H_ot(c,1) =  -mean(log(P_ot_c_cond_past_temp(c,min_order:end)),2,'omitnan');
end





for c = 1:C
    Y_in_temp = Y_in;
    Y_in_temp{c,1} = nan*Y_in_temp{c,1};
    [P_ot_c_cond_past_temp ,alpha_out , b_c_ot_nc_out , P_observ_cond_to_state_out , P_observ_cond_to_state_comp_out]...
        = mc_lsim( pi_steady_lsim , coupling_tetha_IM , transition_matrices_IM ,  gmm_para_lsim  ,  Y_in_temp  );
    
%     [P_ot_c_cond_past_temp ,alpha_out , b_c_ot_nc_out , P_observ_cond_to_state_out , P_observ_cond_to_state_comp_out]...
%         = mc_partial_lsim( pi_steady_lsim , coupling_tetha_IM , transition_matrices_IM ,  gmm_para_lsim  ,  Y_in_temp , c  );
    H_ot_Opast_partial(c,:) =  (-mean(log(P_ot_c_cond_past_temp(:,min_order:end)),2,'omitnan'))';
    
end




PTE = H_ot_Opast_partial - repmat(H_ot_Opast',C,1);
PTE = PTE - diag(diag(PTE));
I = H_ot - H_ot_Opast;
PMI = PTE./repmat(I',C,1);
sig_te = exp(log(sqrt(2*pi*exp(1))) - PTE)/sqrt(2*pi*exp(1));

NCS = PMI.*(1-sig_te);

end

function [P_ot_c_cond_past ,alpha_out , b_c_ot_nc_out , P_observ_cond_to_state_out , P_observ_cond_to_state_comp_out]...
    = mc_lsim( pi_0_lsim , coupling_tetha_IM , transition_matrices_IM ,  gmm_para_lsim  ,  Y_in  )

% Y is observations which is a cell size Cx1  where Y{c,1} is the c'th subsystem observations

% C is the total number of subsystems

%Transition probabilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transitions_matrices is a cell array CxC which cell (i,j) contain
% parameters P(v^j_t|v^i_t-1) (similar to Rabiner) (i->j)

%each cell (i,j) contain a matrix M(i)xM(j) which shows the transition
%probabilities subsystem i on subsystem j according to Rabiner notations  A=Transitions_matrices{i,j} is transition matrix A(n_j , n_i)= P(v^j_t = n_j | v^i_t = n_i ) & sum_n_j (A(i , j)) = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Coupling probabilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coupling matrice is a CxC  which  (c,zee) contain
% parameters Coupling_Tetha(c,zee)
%So sum(Coupling_Tetha)=ones(1,C)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Initials hidden states probabilities of CHMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PI_0{j} contain a vector with length M(j) which contain initials hidden
% states probabilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%emmision probabilities of CHMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHMM_GMM_Param{j}.gmm_para(n_j).P is  mixture wieght (probability) of the
% k'th GMM component related to hidden state (v^j_t = n_j)

% CHMM_GMM_Param{j}.gmm_para(n_j).mu(k).x is  mean vector of the k'th GMM component related to hidden state (v^j_t = n_j)

% CHMM_GMM_Param{j}.gmm_para(n_j).sigma(k).x is  Covariance matrix of the k'th GMM component related to hidden state (v^j_t = n_j)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initializing all parameters
C = size(Y_in , 1);

dim_observation = zeros(C , 1);
channel_num_states = zeros(C , 1);
num_gmm_component =  zeros(C , 1);

for zee = 1:C
    
    dim_observation(zee,1)  = size( Y_in{zee,1} ,1);
    channel_num_states(zee,1)  = size( pi_0_lsim{zee,1} ,1);
    num_gmm_component(zee,1)  = length(gmm_para_lsim{zee,1}.gmm_para(1).P);
    
end

state_numbers_index =  ( [0;cumsum(channel_num_states(:))] );
dimension_numbers_index =  ( [0;cumsum(dim_observation(:))] );

% MU_all is N*C x Dimension
P_all = zeros( C , max(channel_num_states) , max(num_gmm_component)  );
mu_all = zeros(  sum(dim_observation)  , max(channel_num_states) , max(num_gmm_component)   );
sigma_all =  ones(  sum(dim_observation)  , max(channel_num_states) , max(num_gmm_component) );

transitions_matrices = ones( sum(channel_num_states) , sum(channel_num_states)  );
pi_0 = zeros( sum(channel_num_states) , 1  );

for zee = 1:C
    
    gmm_para = gmm_para_lsim{zee,1}.gmm_para ;
    
    for i=1:channel_num_states(zee)
        
        for k=1:num_gmm_component(zee)
            
            P_all(zee,i,k) =  gmm_para(i).P(k,1);
            sigma_all( dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1) ,  i , k) = gmm_para(i).sigma(k).x ;
            mu_all( dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1) ,  i , k) =  gmm_para(i).mu(k).x;
            
        end
        
    end
    
    temp_PI_0 = 0.001+pi_0_lsim{zee,1};
    temp_PI_0 = temp_PI_0/sum(temp_PI_0);
    % initial probabilities
    pi_0(state_numbers_index(zee)+1:state_numbers_index(zee+1) , 1) = temp_PI_0;
    
    for c=1:C
        % Transition probabilities  M(c) x M(zee)
        transitions_matrices( state_numbers_index(c)+1:state_numbers_index(c+1) , state_numbers_index(zee)+1:state_numbers_index(zee+1) ) = transition_matrices_IM{c,zee}  ;
    end
    
end


%Computing Alpha^c_(t|t-1) for all subsystems

channel_time_series = zeros(sum(dim_observation) , size(Y_in{1,1},2));
%Y is PxT, P is observation dimension & T is observations length
for zee=1:C
    channel_time_series(dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1),:)= Y_in{zee,1};
end

T =  ( size( channel_time_series , 2) );

P_ot_c_cond_past = zeros( C , T  );
b_c_ot_nc  = zeros( state_numbers_index(end) , T );

alpha = zeros( state_numbers_index(end) , T  );

%%%beta2= zeros( state_numbers_index(end) , T );
a_alpha_b = zeros(state_numbers_index(end) , C , T-1 );


[P_observ_cond_to_state  ,   P_observ_cond_to_state_comp ]  = ...
    gmm_pdf_fast( P_all , mu_all  , sigma_all  , channel_time_series   ,  channel_num_states ,  dimension_numbers_index , num_gmm_component  );


% zee_index{zee,1} = cell{C,1};
mat_mult = zeros( state_numbers_index(end) , C  );
cuopling_repmat = zeros( state_numbers_index(end) , C  );
cuopling_repmat_beta = zeros( state_numbers_index(end) , C  );
weights_subsystems = zeros( state_numbers_index(end) , C  );
weights_subsystems_opt = weights_subsystems;

% P_vz_tm1_vw_t_O = zeros(state_numbers_index(end) , state_numbers_index(end)  , T );
% P_vz_tm1_vw_t_O_cond = zeros(state_numbers_index(end) , state_numbers_index(end)  , T );

coupling_transmat = zeros( state_numbers_index(end) , state_numbers_index(end)  );
diag_zero = ones( state_numbers_index(end) , state_numbers_index(end)  );


mat_mult_one = mat_mult;
vec_mat_mult = [];

for zee =  1:C
    
    zee_index = ( state_numbers_index(zee)+1:state_numbers_index(zee+1) );
    
    vec_mat_mult = cat(1,vec_mat_mult, ((zee-1)*state_numbers_index(end) + ( state_numbers_index(zee)+1:state_numbers_index(zee+1) ) )' );
    
    zee_coupling = coupling_tetha_IM(:,zee)';
    cuopling_repmat( zee_index ,: ) = repmat( zee_coupling ,channel_num_states(zee),1);
    
    diag_zero( zee_index ,  zee_index )= 0 ;
    
    zee_coupling = coupling_tetha_IM(zee,:);
    cuopling_repmat_beta( zee_index ,: ) = repmat( zee_coupling ,channel_num_states(zee),1);
    
    zee_coupling = coupling_tetha_IM(zee,:);%/sum(coupling_tetha_convex_comb(zee,:));
    weights_subsystems_opt( zee_index ,: ) = repmat( zee_coupling ,channel_num_states(zee),1);
    
    
end


%%%weights_subsystems_opt(isnan(weights_subsystems_opt)) = 1/C;


for zee = 1:C
    zee_index = ( state_numbers_index(zee)+1:state_numbers_index(zee+1) );
    coupling_transmat( : ,  zee_index ) = repmat(cuopling_repmat_beta(:,zee) , 1 , length(zee_index) ).* transitions_matrices(: , zee_index);
end


mat_mult_one(vec_mat_mult)=1;
alpha( :  , 1) = pi_0(:);

P_ot_c_cond_past(: , 1) = mat_mult_one'*(alpha(: , 1).* P_observ_cond_to_state(:,1));
b_c_ot_nc( : , 1) = P_observ_cond_to_state(:,1) ./ (mat_mult_one*P_ot_c_cond_past(: , 1));

%Compute Forward variables & Scaling coefficients
for t = (2:T)
    
    alpha_tm1_tm1 = alpha(: , t-1).*b_c_ot_nc(:,t-1);
    mat_mult(vec_mat_mult) = alpha_tm1_tm1;
    
    a_alpha_b( : , : ,t-1)= transitions_matrices'*mat_mult ;%compute from t=1 to t=T-1
    a_alpha_b_coupled = a_alpha_b(  : , : ,t-1).*cuopling_repmat;
    alpha(: , t) = sum( a_alpha_b_coupled , 2);
    
    P_ot_c_cond_past(: , t) = mat_mult_one'*(alpha(: , t).* P_observ_cond_to_state(:,t));
    
    emmbed_P_ot_c_cond_past = (mat_mult_one*P_ot_c_cond_past(: , t));
    b_c_ot_nc( : , t) = P_observ_cond_to_state(:,t) ./ emmbed_P_ot_c_cond_past;
    
end



if(nargout>1)
    
    alpha_out = cell(C,1);
    b_c_ot_nc_out = cell(C,1);
    P_observ_cond_to_state_out = cell(C,1);
    P_observ_cond_to_state_comp_out = cell(C,1);
    
    for zee = 1:C
        
        alpha_out{zee,1} =  alpha(state_numbers_index(zee)+1:state_numbers_index(zee+1) , :);
        b_c_ot_nc_out{zee,1} =  b_c_ot_nc(state_numbers_index(zee)+1:state_numbers_index(zee+1) , :);
        P_observ_cond_to_state_out{zee,1} =  P_observ_cond_to_state(state_numbers_index(zee)+1:state_numbers_index(zee+1) , :);
        P_observ_cond_to_state_comp_out{zee,1} = P_observ_cond_to_state_comp(state_numbers_index(zee)+1:state_numbers_index(zee+1) , 1:num_gmm_component(zee) , :);
        
    end
    
    
end

end

function [P_ot_c_cond_past ,alpha_out , b_c_ot_nc_out , P_observ_cond_to_state_out , P_observ_cond_to_state_comp_out]...
    = mc_partial_lsim( pi_steady_lsim , coupling_tetha_IM , transition_matrices_IM ,  gmm_para_lsim  ,  Y_in ,zee_in )

% Y is observations which is a cell size Cx1  where Y{c,1} is the c'th subsystem observations

% C is the total number of subsystems

%Transition probabilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transitions_matrices is a cell array CxC which cell (i,j) contain
% parameters P(v^j_t|v^i_t-1) (similar to Rabiner) (i->j)

%each cell (i,j) contain a matrix M(i)xM(j) which shows the transition
%probabilities subsystem i on subsystem j according to Rabiner notations  A=Transitions_matrices{i,j} is transition matrix A(n_j , n_i)= P(v^j_t = n_j | v^i_t = n_i ) & sum_n_j (A(i , j)) = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Coupling probabilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coupling matrice is a CxC  which  (c,zee) contain
% parameters Coupling_Tetha(c,zee)
%So sum(Coupling_Tetha)=ones(1,C)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Initials hidden states probabilities of CHMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PI_0{j} contain a vector with length M(j) which contain initials hidden
% states probabilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%emmision probabilities of CHMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHMM_GMM_Param{j}.gmm_para(n_j).P is  mixture wieght (probability) of the
% k'th GMM component related to hidden state (v^j_t = n_j)

% CHMM_GMM_Param{j}.gmm_para(n_j).mu(k).x is  mean vector of the k'th GMM component related to hidden state (v^j_t = n_j)

% CHMM_GMM_Param{j}.gmm_para(n_j).sigma(k).x is  Covariance matrix of the k'th GMM component related to hidden state (v^j_t = n_j)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initializing all parameters
C = size(Y_in , 1);

dim_observation = zeros(C , 1);
channel_num_states = zeros(C , 1);
num_gmm_component =  zeros(C , 1);

for zee = 1:C
    
    dim_observation(zee,1)  = size( Y_in{zee,1} ,1);
    channel_num_states(zee,1)  = size( pi_steady_lsim{zee,1} ,1);
    num_gmm_component(zee,1)  = length(gmm_para_lsim{zee,1}.gmm_para(1).P);
    
end

state_numbers_index =  ( [0;cumsum(channel_num_states(:))] );
dimension_numbers_index =  ( [0;cumsum(dim_observation(:))] );

% MU_all is N*C x Dimension
P_all = zeros( C , max(channel_num_states) , max(num_gmm_component)  );
mu_all = zeros(  sum(dim_observation)  , max(channel_num_states) , max(num_gmm_component)   );
sigma_all =  ones(  sum(dim_observation)  , max(channel_num_states) , max(num_gmm_component) );

transitions_matrices = ones( sum(channel_num_states) , sum(channel_num_states)  );
pi_0 = zeros( sum(channel_num_states) , 1  );

for zee = 1:C
    
    gmm_para = gmm_para_lsim{zee,1}.gmm_para ;
    
    for i=1:channel_num_states(zee)
        
        for k=1:num_gmm_component(zee)
            
            P_all(zee,i,k) =  gmm_para(i).P(k,1);
            sigma_all( dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1) ,  i , k) = gmm_para(i).sigma(k).x ;
            mu_all( dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1) ,  i , k) =  gmm_para(i).mu(k).x;
            
        end
        
    end
    
    temp_PI_0 = pi_steady_lsim{zee,1};
    temp_PI_0 = temp_PI_0/sum(temp_PI_0);
    % initial probabilities
    pi_0(state_numbers_index(zee)+1:state_numbers_index(zee+1) , 1) = temp_PI_0;
    
    for c=1:C
        % Transition probabilities  M(c) x M(zee)
        transitions_matrices( state_numbers_index(c)+1:state_numbers_index(c+1) , state_numbers_index(zee)+1:state_numbers_index(zee+1) ) = transition_matrices_IM{c,zee}  ;
    end
    
end


%Computing Alpha^c_(t|t-1) for all subsystems

channel_time_series = zeros(sum(dim_observation) , size(Y_in{1,1},2));
%Y is PxT, P is observation dimension & T is observations length
for zee=1:C
    channel_time_series(dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1),:)= Y_in{zee,1};
end

T =  ( size( channel_time_series , 2) );

P_ot_c_cond_past = zeros( C , T  );
b_c_ot_nc  = zeros( state_numbers_index(end) , T );

alpha = zeros( state_numbers_index(end) , T  );

%%%beta2= zeros( state_numbers_index(end) , T );
a_alpha_b = zeros(state_numbers_index(end) , C , T-1 );


[P_observ_cond_to_state  ,   P_observ_cond_to_state_comp ]  = ...
    gmm_pdf_fast( P_all , mu_all  , sigma_all  , channel_time_series   ,  channel_num_states ,  dimension_numbers_index , num_gmm_component  );


% zee_index{zee,1} = cell{C,1};
mat_mult = zeros( state_numbers_index(end) , C  );
cuopling_repmat = zeros( state_numbers_index(end) , C  );
cuopling_repmat_beta = zeros( state_numbers_index(end) , C  );
weights_subsystems = zeros( state_numbers_index(end) , C  );
weights_subsystems_opt = weights_subsystems;

% P_vz_tm1_vw_t_O = zeros(state_numbers_index(end) , state_numbers_index(end)  , T );
% P_vz_tm1_vw_t_O_cond = zeros(state_numbers_index(end) , state_numbers_index(end)  , T );

coupling_transmat = zeros( state_numbers_index(end) , state_numbers_index(end)  );
diag_zero = ones( state_numbers_index(end) , state_numbers_index(end)  );


mat_mult_one = mat_mult;
vec_mat_mult = [];

for zee =  1:C
    
    zee_index = ( state_numbers_index(zee)+1:state_numbers_index(zee+1) );
    
    vec_mat_mult = cat(1,vec_mat_mult, ((zee-1)*state_numbers_index(end) + ( state_numbers_index(zee)+1:state_numbers_index(zee+1) ) )' );
    
    zee_coupling = coupling_tetha_IM(:,zee)';
    cuopling_repmat( zee_index ,: ) = repmat( zee_coupling ,channel_num_states(zee),1);
    
    diag_zero( zee_index ,  zee_index )= 0 ;
    
    zee_coupling = coupling_tetha_IM(zee,:);
    cuopling_repmat_beta( zee_index ,: ) = repmat( zee_coupling ,channel_num_states(zee),1);
    
    zee_coupling = coupling_tetha_IM(zee,:);%/sum(coupling_tetha_convex_comb(zee,:));
    weights_subsystems_opt( zee_index ,: ) = repmat( zee_coupling ,channel_num_states(zee),1);
    
    
end


%%%weights_subsystems_opt(isnan(weights_subsystems_opt)) = 1/C;


for zee = 1:C
    zee_index = ( state_numbers_index(zee)+1:state_numbers_index(zee+1) );
    coupling_transmat( : ,  zee_index ) = repmat(cuopling_repmat_beta(:,zee) , 1 , length(zee_index) ).* transitions_matrices(: , zee_index);
end


mat_mult_one(vec_mat_mult)=1;
alpha( :  , 1) = pi_0(:);

P_ot_c_cond_past(: , 1) = mat_mult_one'*(alpha(: , 1).* P_observ_cond_to_state(:,1));
b_c_ot_nc( : , 1) = P_observ_cond_to_state(:,1) ./ (mat_mult_one*P_ot_c_cond_past(: , 1));

%Compute Forward variables & Scaling coefficients
for t = (2:T)
    
    alpha_tm1_tm1 = alpha(: , t-1).*b_c_ot_nc(:,t-1);
    mat_mult(vec_mat_mult) = alpha_tm1_tm1;
    
    a_alpha_b( : , : ,t-1)= transitions_matrices'*mat_mult ;%compute from t=1 to t=T-1
    a_alpha_b_coupled = a_alpha_b(  : , : ,t-1).*cuopling_repmat;
    alpha(: , t) = sum( a_alpha_b_coupled , 2);
    alpha(state_numbers_index(zee_in)+1:state_numbers_index(zee_in+1) , t) = pi_steady_lsim{zee_in,1};
    P_ot_c_cond_past(: , t) = mat_mult_one'*(alpha(: , t).* P_observ_cond_to_state(:,t));
    
    emmbed_P_ot_c_cond_past = (mat_mult_one*P_ot_c_cond_past(: , t));
    b_c_ot_nc( : , t) = P_observ_cond_to_state(:,t) ./ emmbed_P_ot_c_cond_past;
    
end



if(nargout>1)
    
    alpha_out = cell(C,1);
    b_c_ot_nc_out = cell(C,1);
    P_observ_cond_to_state_out = cell(C,1);
    P_observ_cond_to_state_comp_out = cell(C,1);
    
    for zee = 1:C
        
        alpha_out{zee,1} =  alpha(state_numbers_index(zee)+1:state_numbers_index(zee+1) , :);
        b_c_ot_nc_out{zee,1} =  b_c_ot_nc(state_numbers_index(zee)+1:state_numbers_index(zee+1) , :);
        P_observ_cond_to_state_out{zee,1} =  P_observ_cond_to_state(state_numbers_index(zee)+1:state_numbers_index(zee+1) , :);
        P_observ_cond_to_state_comp_out{zee,1} = P_observ_cond_to_state_comp(state_numbers_index(zee)+1:state_numbers_index(zee+1) , 1:num_gmm_component(zee) , :);
        
    end
    
    
end

end


