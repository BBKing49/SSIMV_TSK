function [best_pg, best_v, best_b, best_TSK_acc_te, best_TSK_acc_std] = train_TSK_FS( X_t_view , labels, folds_num)
%X_t_view Ϊt�ӽ����� N*D��YΪ��ǩ N*1, ��C�࣬��TΪN*C,pgΪM*C,MΪ������,�ⲿ��һ��
Ms = 2:2:10;
lamdas = [2.^(-5:5)];
T = lab2vec(labels);
clusters = labels;
masks_te=cv_fold(folds_num,T);
best_TSK_acc_te = 0;
for lamda = lamdas
    for M = Ms
        %5��������֤Ѱ��
        result = zeros(folds_num,1);
        for fold=1:folds_num
            X = X_t_view;
            mask_te=masks_te{fold,1};
            mask_tr=~mask_te;
            X_tr=X(mask_tr,:);
            T_tr=T(mask_tr,:);
            X_te=X(mask_te,:);
            clusters_te = clusters(mask_te,:);
            %ѵ��
            [v,b] = gene_ante_deter(X_tr,M);
            Xg = fromXtoZ(X_tr,v,b);   %Xg:N*K
            Xg1 = Xg'*Xg;
            pg = pinv(Xg1 + lamda*eye( size(Xg1)))*Xg'*T_tr;   
            %����
            [Y_te] = test_TSK_FS( X_te , pg, v, b);
            labels_te = vec2lab(Y_te);
            acc_mean=sum(labels_te==clusters_te)/length(clusters_te);
            result(fold,1)=acc_mean;
        end
        acc_te_mean = mean(result(:,1));
        acc_te_std = std(result(:,1));
        if acc_te_mean>best_TSK_acc_te
            best_TSK_acc_te = acc_te_mean;
            best_TSK_acc_std = acc_te_std;
            best_pg = pg;
            best_v = v;
            best_b = b;
        end
    end %end M
end %end lamda
