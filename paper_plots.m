clear auroc_lsim auroc_hermes
clc
close all


% T_durations = [100,200,500,1000];
% db_vals = -10:5:20;


load('results_fast2.mat')

save_path = './figs';

if ~exist(save_path, 'dir')
    mkdir(save_path);
end

gaussian_noise = 1;

pte_lsim = results.pte_lsim;
pte_hermes = results.pte_hermes;

system_types = 1;

for system_type = system_types

    for Npop = channels(system_type).num
        count_T = 0;
        for T = T_durations
            count_T = count_T + 1;
            db_count = 0;
            for  db_quality = db_vals
                close all
                db_count = db_count + 1;

                G = G_all{system_type,count_T,Npop,db_count,gaussian_noise};
                est_G = est_G_all{system_type,count_T,Npop,db_count,gaussian_noise} ;
                est_G2 = est_G2_all{system_type,count_T,Npop,db_count,gaussian_noise} ;
                G = double(G);
                for n1=1:Npop
                    G(n1,n1,:)=nan;
                end

                sen(Npop) = sum((est_G(:)==G(:))&(G(:)==1),'omitnan')/sum(G(:)==1);% recall
                precision(Npop) = sum((est_G(:)==G(:))&(est_G(:)==1),'omitnan')/sum(est_G(:)==1);
                spec(Npop) = sum((est_G(:)==G(:))&(G(:)==0),'omitnan')/sum(G(:)==0);
                Acc(system_type,count_T,Npop,db_count,gaussian_noise) = sum((est_G(:)==G(:)),'omitnan')/sum(~isnan(G(:)));
                fscore(system_type,count_T,Npop,db_count,gaussian_noise) = 2*(sen(Npop)*precision(Npop))/(10^-5+sen(Npop)+precision(Npop));
                bacc(system_type,count_T,Npop,db_count,gaussian_noise) = (sen(Npop) +  spec(Npop) )/2;

                sen2(Npop) = sum((est_G2(:)==G(:))&(G(:)==1),'omitnan')/sum(G(:)==1);
                precision2(Npop) = sum((est_G2(:)==G(:))&(est_G2(:)==1),'omitnan')/sum(est_G2(:)==1);
                spec2(Npop) = sum((est_G2(:)==G(:))&(G(:)==0),'omitnan')/sum(G(:)==0);
                Acc2(system_type,count_T,Npop,db_count,gaussian_noise) = sum((est_G2(:)==G(:)),'omitnan')/sum(~isnan(G(:)));
                fscore2(system_type,count_T,Npop,db_count,gaussian_noise) = 2*(sen2(Npop)*precision2(Npop))/(10^-5+sen2(Npop)+precision2(Npop));
                bacc2(system_type,count_T,Npop,db_count,gaussian_noise) = (sen2(Npop) +  spec2(Npop) )/2;


                G = G_all{system_type,count_T,Npop,db_count,gaussian_noise};
                G = double(G);
                for n1=1:size(G,1)
                    G(n1,n1,:)=nan;
                end


                PTE_org = pte_lsim{system_type,count_T,Npop,db_count,gaussian_noise} ;
                PTE_hermes = pte_hermes{system_type,count_T,Npop,db_count,gaussian_noise};


                est_lsim = PTE_org(:);
                est_hermes = PTE_hermes(:);

                G_true = G(:);
                indnan = isnan(G_true);
                est_lsim(indnan)=[];
                est_hermes(indnan)=[];
                %                 est_lsim(est_lsim<0)=0;
                %                 est_hermes(est_hermes<0)=0;
                G_true(indnan)=[];

                try
                    %[X,Y,TT,auroc_lsim(system_type,count_T,Npop,db_count,gaussian_noise)] = perfcurve(G_true,est_lsim,1);
                    %[X2,Y2,TT2,auroc_hermes(system_type,count_T,Npop,db_count,gaussian_noise)] = perfcurve(G_true,est_hermes,1);
                    %[prec, tpr, fpr, thresh, auroc_lsim(system_type,count_T,Npop,db_count,gaussian_noise),aucpr_lsim(system_type,count_T,Npop,db_count,gaussian_noise)] = prec_rec(est_lsim,G_true,'plotPR',1,'plotROC',1,'numThresh',length(G_true));
                    %[prec, tpr, fpr, thresh, auroc_hermes(system_type,count_T,Npop,db_count,gaussian_noise),aucprc_hermes(system_type,count_T,Npop,db_count,gaussian_noise)] = prec_rec(est_hermes,G_true,'plotPR',1,'plotROC',1,'numThresh',length(G_true),'holdFigure',1);

                    [prec2, tpr2, fpr2, thresh, auroc_hermes(system_type,count_T,Npop,db_count,gaussian_noise),aucprc_hermes(system_type,count_T,Npop,db_count,gaussian_noise)] = prec_rec(est_hermes,G_true,'numThresh',min(1000,length(G_true)));
                    [prec1, tpr1, fpr1, thresh, auroc_lsim(system_type,count_T,Npop,db_count,gaussian_noise),aucpr_lsim(system_type,count_T,Npop,db_count,gaussian_noise)] = prec_rec(est_lsim,G_true,'numThresh',min(1000,length(G_true)));
                catch
                    auroc_lsim(system_type,count_T,Npop,db_count,gaussian_noise) = nan;
                    auroc_hermes(system_type,count_T,Npop,db_count,gaussian_noise) = nan;

                    aucpr_lsim(system_type,count_T,Npop,db_count,gaussian_noise) = nan;
                    aucprc_hermes(system_type,count_T,Npop,db_count,gaussian_noise) = nan;
                end

                if system_type==1 && count_T>=3 && db_count==3
                    close all
                    a = figure;
                    plot(fpr1, tpr1,'b','LineWidth',2)
                    hold on
                    grid on
                    plot(fpr2, tpr2,'r','LineWidth',2)
                    legend({'$LSIM$','$kNN$'},'Location','southeast','FontSize',15,'Interpreter' ,'latex')
                    set(gca, 'FontWeight','bold','FontSize',14);
                    xlabel('FPR','FontSize',17,'Interpreter' ,'latex')
                    ylabel('TPR','FontSize',17,'Interpreter' ,'latex')
                    ylim([0,1])
                    title(['$a.~ROC~curve$'],'FontSize',17,'Interpreter' ,'latex')
                    saveas(a,[save_path,'/dauroc_sys',num2str(system_type),'_T',num2str(count_T),'_C',num2str(Npop),'_G',num2str(gaussian_noise),'.eps'],'epsc2')

                    a = figure;
                    plot(tpr1, prec1,'b','LineWidth',2)
                    hold on
                    grid on
                    plot(tpr2, prec2,'r','LineWidth',2)
                    legend({'$LSIM$','$kNN$'},'Location','southeast','FontSize',15,'Interpreter' ,'latex')
                    set(gca, 'FontWeight','bold','FontSize',14);
                    xlabel('Recall','FontSize',17,'Interpreter' ,'latex')
                    ylabel('Precision','FontSize',17,'Interpreter' ,'latex')
                    ylim([0,1])
                    title(['$PR~curve$'],'FontSize',17,'Interpreter' ,'latex')
                    saveas(a,[save_path,'/dauprc_sys',num2str(system_type),'_T',num2str(count_T),'_C',num2str(Npop),'_G',num2str(gaussian_noise),'.eps'],'epsc2')
                    saveas(a,[save_path,'/dauprc_sys',num2str(system_type),'_T',num2str(count_T),'_C',num2str(Npop),'_G',num2str(gaussian_noise),'.fig'])

                    yu=0;
                    close all
                end

            end
        end
    end


    bacc(:,:,1,:,:) = nan;
    bacc2(:,:,1,:,:) = nan;
    fscore(:,:,1,:,:) = nan;
    fscore2(:,:,1,:,:) = nan;

    %%
    clc

    %T_durations = [100,200,500,1000];
    % db_vals = -10:5:20

    c=0;
    for channels_num = [3,5,10]
        c=c+1;
        for count_T = 3:4

            close all
            a = figure;
            plot(db_vals, squeeze(auroc_lsim(system_type,count_T,channels_num,:,gaussian_noise)),'bo--','LineWidth',1.5)
            hold on
            grid on
            plot(db_vals, squeeze(auroc_hermes(system_type,count_T,channels_num,:,gaussian_noise)),'rx-.','LineWidth',1.5)
            set(gca, 'FontWeight','bold','FontSize',14);
            xlabel('SNR (dB)','FontSize',17,'Interpreter' ,'latex')
            ylabel('AUC','FontSize',17,'Interpreter' ,'latex')
            legend({'$LSIM$','$kNN$'},'Location','southeast','FontSize',15,'Interpreter' ,'latex')
            if count_T==3
                title(['$a.~T=',num2str(T_durations(count_T)),'$'],'FontSize',17,'Interpreter' ,'latex')
            else
                title(['$b.~T=',num2str(T_durations(count_T)),'$'],'FontSize',17,'Interpreter' ,'latex')
            end
            ylim([0.4,1])
            saveas(a,[save_path,'/auc_sys',num2str(system_type),'_T',num2str(count_T),'_C',num2str(channels_num),'_G',num2str(gaussian_noise),'.jpg'])
            saveas(a,[save_path,'/auc_sys',num2str(system_type),'_T',num2str(count_T),'_C',num2str(channels_num),'_G',num2str(gaussian_noise),'.eps'],'epsc2')

            if count_T==3
                diff_auc_500(c,:) = squeeze(auroc_lsim(system_type,count_T,channels_num,:,gaussian_noise))-squeeze(auroc_hermes(system_type,count_T,channels_num,:,gaussian_noise));
            else
                diff_auc_1000(c,:) = squeeze(auroc_lsim(system_type,count_T,channels_num,:,gaussian_noise))-squeeze(auroc_hermes(system_type,count_T,channels_num,:,gaussian_noise));
            end
            a = figure;
            plot(db_vals, squeeze(aucpr_lsim(system_type,count_T,channels_num,:,gaussian_noise)),'bo--','LineWidth',1.5)
            hold on
            grid on
            plot(db_vals, squeeze(aucprc_hermes(system_type,count_T,channels_num,:,gaussian_noise)),'rx-.','LineWidth',1.5)
            legend({'$LSIM$','$kNN$'},'Location','southeast','FontSize',15,'Interpreter' ,'latex')
            set(gca, 'FontWeight','bold','FontSize',14);
            xlabel('SNR (dB)','FontSize',17,'Interpreter' ,'latex')
            ylabel('AUCPR','FontSize',17,'Interpreter' ,'latex')
            ylim([0,1])
            if count_T==3
                title(['$a.~T=',num2str(T_durations(count_T)),'$'],'FontSize',17,'Interpreter' ,'latex')
            else
                title(['$b.~T=',num2str(T_durations(count_T)),'$'],'FontSize',17,'Interpreter' ,'latex')
            end
            saveas(a,[save_path,'/auprc_sys',num2str(system_type),'_T',num2str(count_T),'_C',num2str(channels_num),'_G',num2str(gaussian_noise),'.fig'])
            saveas(a,[save_path,'/auprc_sys',num2str(system_type),'_T',num2str(count_T),'_C',num2str(channels_num),'_G',num2str(gaussian_noise),'.eps'],'epsc2')
            if count_T==3
                diff_aucpr_500(c,:) = squeeze(aucpr_lsim(system_type,count_T,channels_num,:,gaussian_noise))-squeeze(aucprc_hermes(system_type,count_T,channels_num,:,gaussian_noise));
            else
                diff_aucpr_1000(c,:) = squeeze(aucpr_lsim(system_type,count_T,channels_num,:,gaussian_noise))-squeeze(aucprc_hermes(system_type,count_T,channels_num,:,gaussian_noise));
            end

        end
    end

end




