function [best_result,  TSK_result, H_cells ] = expt_SSIMV_TSK( data, labels ,incomplete_rate, labeled_rate, output_file)
%数据参数初始化
view_nums = size(data,1); %未添加隐空间的视角数

TSK_cell = cell(view_nums,4); %第一列为pg，第二列为v，第三列为b，第四列为w(权重)
folds_num = 5;
options.view_nums = view_nums;
maxIter = 5;
options.maxIter = maxIter;

lamda1s = 2.^(-5:5);
lamda2s = 2.^(-5:5);
lamda3s = 2.^(-5:5);
lamda4s = 2.^(-5:5);
lamda5s = 2.^(-5:5);

N = numel(labels);
randindex = randperm(N);
labels = labels(randindex,:);
T = lab2vec(labels);
c = size(T,2);

for i=1:view_nums
    data{i} = data{i}(randindex,:);
end

%随机删去一部分数据，生成不完整多视角数据
p = randperm(numel(labels)-1);
p = p+1;
p = repmat(p, 1,view_nums);
incomplete_size = ceil(incomplete_rate*numel(labels)*(1/view_nums));

for i = 1:view_nums
    incomplete_index{i} = p((i-1)*incomplete_size+1:i*incomplete_size);
    incomplete_index{i} = sort(incomplete_index{i});
    tmp = 1:N;tmp=tmp';
    tmp(incomplete_index{i},:) = [];
    complete_index{i} = tmp;
end

for i = 1:view_nums
    
    for j = 1:numel(incomplete_index{i})
        total = sum(data{i}(1:incomplete_index{i}(j)-1,:),1);
        incomplete = sum(data{i}(incomplete_index{i}(1:j-1),:),1);
        data{i,1}(incomplete_index{i}(j),:) = (total - incomplete)./(incomplete_index{i}(j)-j);%0;
    end
    data{i} = mapminmax(data{i}', 0, 1)';    
    W{i} = diag(sparse(ones(N,1)));
    E{i} = diag(zeros(N,1));
    counter = 0;
    for j = 1:numel(incomplete_index{i})
        counter = counter +1;
        E{i}(incomplete_index{i}(j),incomplete_index{i}(j)) = 1;
        W{i}(incomplete_index{i}(j),incomplete_index{i}(j)) = 1.0*counter/incomplete_index{i}(j);
    end

end

labeled_index = [];
unlabeled_index = [];
% 划分标记和未标记
for i = 1:c
    ind = find(labels==i);
    tmp = round(length(ind)*labeled_rate);
    labeled_index = [labeled_index;ind(1:tmp)];
    unlabeled_index = [unlabeled_index;ind(tmp+1:end)];
end
labeled_index = sort(labeled_index);
unlabeled_index = sort(unlabeled_index);

%对每个视角训练一个TSK
tic;
TSK_result = zeros(view_nums,2);

for view_num = 1:view_nums
    acc_data = data{view_num,1}(labeled_index,:);
    labels_tmp = labels(labeled_index,:);
    [pg, v, b, single_best_acc, single_best_acc_std] = train_TSK_FS( acc_data , labels_tmp, folds_num);
    TSK_cell{ view_num, 1 } = pg;
    TSK_cell{ view_num, 2 } = v;
    TSK_cell{ view_num, 3 } = b;
    TSK_cell{ view_num, 4 } = 1/view_nums;
    TSK_result(view_num,1) = single_best_acc;
    TSK_result(view_num,2) = single_best_acc_std;
end
toc;

T_tr=T;
T_tr(unlabeled_index,:) = rand(length(unlabeled_index),c);
temp(:,:)=sum(T_tr(unlabeled_index,:),2);
T_tr(unlabeled_index,:) = T_tr(unlabeled_index,:)./temp;

best_acc_te = 0;

for lamda1 = lamda1s
    for lamda2 =  lamda2s
        for lamda3 =  lamda3s
            for lamda4 = lamda4s
                for lamda5 = lamda5s
                    options.lamda1 = lamda1;
                    options.lamda2 = lamda2;
                    options.lamda3 = lamda3;
                    options.lamda4 = lamda4;
                    options.lamda5 = lamda5;
                                      
                    clusters_te = vec2lab(T(unlabeled_index,:));
                    options.unlabeled_index = unlabeled_index;
                    options.incomplete_index = incomplete_index;
                    tic;
                    [best_TSK_cell, data_new, T_trans] = SSIMV_TSK( data, TSK_cell, T_tr, W, E, complete_index, options);
                    t=toc;
%                     labels_trans = vec2lab(T_trans(unlabeled_index,:));
%                     acc_trans=sum(labels_trans==clusters_te)/length(clusters_te);
                    
                    [Y_te] = test_mul_TSK( data_new ,best_TSK_cell, view_nums, c, unlabeled_index);
                    labels_te = vec2lab(Y_te);
                    acc_te=sum(labels_te==clusters_te)/length(clusters_te);

                    if acc_te>best_acc_te
                        best_acc_te = acc_te;
                        best_result.result_pg = best_acc_te;
                        best_result.best_model = best_TSK_cell;
                        
                        best_result.lamda1 = lamda1;
                        best_result.lamda2 = lamda2;
                        best_result.lamda3 = lamda3;
                        best_result.lamda4 = lamda4;
                        best_result.lamda5 = lamda5;
                        best_result.time=t;
                        acc_te
                    end

                end
            end
        end 
    end 
end 

