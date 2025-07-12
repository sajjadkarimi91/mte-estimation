close all
clc
clear

addpath(genpath(pwd))


%%

max_itration = 100;

num_surrogate = 20;
surrogate_data = [];
surro_type = 2; % 1:IAAFT 2: Shift
max_repeat = 100;
T_durations = [100,200,500,1000];
db_vals = -10:5:20;
channels(1).num = [3,5,10];%Neural Mass Models
channels(2).num = 2;%nCREANN
channels(3).num = 2:10;%henon_map
channels(4).num =  3;% nlmvar_exp
channels(5).num = [2:10];% high_dim_network

fs = 100; %100 or 200
system_types = [1];
nms_lags = 1;
n_lags_hermes = 2;
n_lags_lsim = 1;

params.max_delay = nms_lags;
params.fs = fs;
nms_selection_rate = 0.2;

for gaussian_noise = 1

    for system_type = system_types

        for Npop = channels(system_type).num
            count_T = 0;
            for T = T_durations
                count_T = count_T + 1;

                for num_rep = 1:max_repeat
                    params.max_delay = nms_lags;

                    %% 2020, Transfer Entropy as a Measure of Brain Connectivity: A Critical Analysis With the Help of Neural Mass Models
                    % 2010,  The generation of rhythms within a cortical region: analysis of a neural mass model

                    % definition of excitatory and inhibitory synapses
                    Wp=zeros(Npop);

                    %excitatorty synapses
                    index_one = rand(Npop,Npop);
                    Wp(index_one>(1-nms_selection_rate/2))= 40;   % Wp_12 = 0, 20, 40, 60, 80
                    Wp = Wp.*(1-eye(Npop));

                    % inhibitory synapses
                    %         Wf = Wp;
                    %         Wf(Wf>0) = 0; % 0, 20, 40, 60 : 4 panels
                    index_one  = rand(Npop,Npop);
                    Wf=zeros(Npop);
                    Wf(index_one>(1-nms_selection_rate/2))=40;   % 0, 20, 40, 60 : 4 panels
                    Wf = Wf.*(1-eye(Npop));

                    params.Wp = Wp;
                    params.Wf = Wf;
                    params.Npop = Npop;


                    [eeg, G_eff, fs, tt] = nms_connectivity_data(T, params);
                    multichannel_timeseries_sim = eeg;



                    multichannel_timeseries_sim_org = multichannel_timeseries_sim;


                    db_count = 0;
                    for  db_quality = db_vals
                        clc
                        system_type
                        Npop
                        T
                        num_rep
                        db_quality
                        db_count = db_count+1;

                        lambda = 10^(-db_quality/20);
                        if gaussian_noise == 1
                            noise1 = randn(size(multichannel_timeseries_sim_org));
                        else
                            noise1 = randn(size(multichannel_timeseries_sim_org))+3;
                            noise2 = randn(size(multichannel_timeseries_sim_org))-3;
                            ind_perm = randperm(length(noise1(:)));
                            noise1(ind_perm(1:length(noise1(:))/2)) = noise2(ind_perm(1:length(noise1(:))/2));
                        end

                        multichannel_timeseries_sim = multichannel_timeseries_sim_org + lambda*norm(multichannel_timeseries_sim_org(:))/norm(noise1(:))* noise1;


                        if Npop==2
                            config.measures={'TE'};
                        else
                            config.measures={'PTE'};
                        end
                        config.window.length = 1000*T/fs; %window: Window length in milliseconds.
                        config.window.overlap = 0; %Percentage
                        config.window.fs = fs;
                        config.TimeDelay = 1;
                        config.EmbDim = 1;
                        config.Nneighbours = ceil(sqrt(T));
                        config.Nlags = n_lags_hermes;
                        try
                        output_IT = H_methods_IT_corrected ( multichannel_timeseries_sim', config );
                        catch ME
                            disp('You need mex-file "tim_matlab" be enabled from Hermes Toolbox');
                        end
                        if Npop==2
                            PTE_hermes = output_IT.TE.data;
                        else
                            PTE_hermes = output_IT.PTE.data;
                        end



                        % lsim modelling
                        clear channels_observations index_channels index_lags
                        counter_c = 0;
                        for c=1:size(multichannel_timeseries_sim,1)
                            for n = 0:n_lags_lsim-1
                                counter_c =counter_c+1;
                                channels_observations{counter_c,1} = multichannel_timeseries_sim(c,1+n_lags_lsim-n:end-n_lags_lsim-n);
                                index_channels(counter_c) = c;
                                index_lags(counter_c) = n;

                            end
                        end

                        clear Log Model_t Model_rep BIC_S channel_num_states num_gmm_component Model_lsim AIC_lsim_S


                        if num_rep==1
                            state_max(db_count) = 8;
                            gmm_max(db_count) = 2;
                            state_min = 2;
                        else
                            state_min = max(2,state_max(db_count)-2);
                        end


                        state_numbers_all = state_min:state_max(db_count);

                        if gmm_max(db_count) > 2
                            num_gmm_component_all = [ones(1,length(state_numbers_all)),2*ones(1,length(state_numbers_all)),3*ones(1,length(state_numbers_all))];
                            state_numbers_all = [state_numbers_all,state_numbers_all,state_numbers_all];
                        elseif gmm_max(db_count) > 1
                            num_gmm_component_all = [ones(1,length(state_numbers_all)),2*ones(1,length(state_numbers_all))];
                            state_numbers_all = [state_numbers_all,state_numbers_all];
                        else
                            num_gmm_component_all = ones(1,length(state_numbers_all));
                        end
                        %  state_numbers_all = [2:6];
                        %  num_gmm_component_all = ones(1,length(state_numbers_all));
                        %  num_gmm_component_all = [ones(1,length(state_numbers_all)),2*ones(1,length(state_numbers_all))];
                        %  state_numbers_all = [state_numbers_all,state_numbers_all];

                        C = size(channels_observations,1);
                        extra.plot = 0;
                        extra.check_convergence=0;

                        parfor s_search = 1: length(state_numbers_all)

                            channel_num_states = ones(1,C)*state_numbers_all(s_search);
                            num_gmm_component = ones(1,C)*num_gmm_component_all(s_search);

                            [~ , coupling_tetha_convex_comb , transition_matrices_convex_comb ,  lsim_gmm_para ,  AIC , log_likelihood , BIC ,pi_steady] = ...
                                em_lsim( channels_observations , channel_num_states , num_gmm_component , max_itration , extra);

                            Model_lsim{s_search}.CHMM_GMM_Param = lsim_gmm_para ;
                            Model_lsim{s_search}.Transitions_matrices = transition_matrices_convex_comb;
                            Model_lsim{s_search}.Coupling_Tetha = coupling_tetha_convex_comb;
                            Model_lsim{s_search}.PI_0 = pi_steady ;
                            BIC_lsim_S(s_search) = BIC(end);
                            AIC_lsim_S(s_search) = AIC(end);

                        end

                        [~,Index_min] = min(AIC_lsim_S);
                        lsim_gmm_para = Model_lsim{Index_min}.CHMM_GMM_Param ;
                        transition_matrices_convex_comb = Model_lsim{Index_min}.Transitions_matrices;
                        coupling_tetha_convex_comb = Model_lsim{Index_min}.Coupling_Tetha;
                        pi_0_lsim = Model_lsim{Index_min}.PI_0 ;

                        state_max(db_count) = length(pi_0_lsim{1}) + 1;
                        gmm_max(db_count) = num_gmm_component_all(Index_min);
                        extra_te.flag_exact = 0;
                        extra_te.min_order = 10;
                        extra_te.index_channels = index_channels;
                        extra_te.index_lags = index_lags;
                        [PTE_org ,I , PMI , NCS , sig_te] = te_embed_lsim( pi_0_lsim , coupling_tetha_convex_comb , transition_matrices_convex_comb ,  lsim_gmm_para  ,  channels_observations, extra_te);

                        channel_num_states = ones(1,C)*state_numbers_all(Index_min);
                        num_gmm_component = ones(1,C)*num_gmm_component_all(Index_min);


                        pte_lsim_surrogate = std(PTE_org(:),'omitnan')*randn(Npop,Npop,num_surrogate);
                        pte_hermes_surro = std(PTE_hermes(:),'omitnan')*randn(Npop,Npop,num_surrogate);


                        pte_hermes_surro(isnan(pte_lsim_surrogate)) = nan;

                        pte_lsim_surrogate(pte_lsim_surrogate<0)=0;
                        mn_PTE1 = mean(pte_lsim_surrogate(:),'omitnan');
                        std_PTE1 = std(pte_lsim_surrogate(:),'omitnan');
                        est_G_lsim = (mn_PTE1 + 3*std_PTE1) < PTE_org;

                        pte_hermes_surro(pte_hermes_surro<0)=0;
                        mn_PTE = mean(pte_hermes_surro(:),'omitnan');
                        std_PTE = std(pte_hermes_surro(:),'omitnan');
                        est_G_hermes = mn_PTE + 3*std_PTE < PTE_hermes;

                        results.pte_lsim{system_type,count_T,Npop,db_count,gaussian_noise}(:,:,num_rep) = PTE_org;
                        results.pte_hermes{system_type,count_T,Npop,db_count,gaussian_noise}(:,:,num_rep) = PTE_hermes;

                        results.pte_lsim_surro{system_type,count_T,Npop,db_count,gaussian_noise}(:,num_rep) = pte_lsim_surrogate(:);
                        results.pte_hermes_surro{system_type,count_T,Npop,db_count,gaussian_noise}(:,num_rep) = pte_hermes_surro(:);

                        %   G_eff
                        %   est_G_lsim
                        %   est_G_hermes
                        G_all{system_type,count_T,Npop,db_count,gaussian_noise}(:,:,num_rep) = G_eff;
                        est_G_all{system_type,count_T,Npop,db_count,gaussian_noise}(:,:,num_rep) = est_G_lsim;
                        est_G2_all{system_type,count_T,Npop,db_count,gaussian_noise}(:,:,num_rep) = est_G_hermes;
                    end
                end
            end
            save('nms_mte_results.mat',"est_G2_all","est_G_all","G_all",'results','T_durations', 'db_vals' , 'system_types' , 'channels')

        end

    end

end