function [v,b] = gene_ante_deter(data,K)
% 2018-05-24 PengXu x1724477385@126.com
% Generate the v and b deterministically.
k = K;
C = var_part(data,k);
v = C;
b = kernel_width(C,data,0,1);

end

function [C] = var_part(data,K)
% Realize the Var-Part algorithm proposed in "In Search of Deterministic 
% Methods for Initializing K-Means and Gaussian Mixture" by Ting Su and
% et.al. 2006

% This algorithm can generate the clustering centers deterministically like
% K-Means but the latter generate the centers randomly.

% data: n_examples * n_features
% K: the number of clusters
% C::return: n_clusters * n_features

clusters = cell(1,K);
C = zeros(K,size(data,2));
k = 1;
while (k<K)
    var_dimen = var(data,0,1);
    [~,maxvar_index] = max(var_dimen);
    data_maxvar = data(:,maxvar_index);
    mean_maxvar = mean(data_maxvar);
    class{1} = data(data_maxvar<=mean_maxvar,:);
    class{2} = data(data_maxvar>mean_maxvar,:);
    sum_norm_1 = scatter_within(class{1});
    sum_norm_2 = scatter_within(class{2});
    [~,max_index_class] = max([sum_norm_1,sum_norm_2]);
    [~,min_index_class] = min([sum_norm_1,sum_norm_2]);
    data = class{max_index_class}; % for partition in next loop
    clusters{1,k} = class{min_index_class};
    k = k+1;
    if k==K
        clusters{1,k} = class{max_index_class};
    end
end

for i=1:K
    C(i,:) = mean(clusters{1,i},1);
end

end

function sum_norm = scatter_within(data)
mean_data = mean(data,1);
mean_data = repmat(mean_data,size(data,1),1);
data_new = data-mean_data;
sum_norm = 0;
for i=1:size(data,1)
    sum_norm = sum_norm+norm(data_new(i,:));
end
end
