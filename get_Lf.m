function [cell_l,cell_s] = get_Lf(X,unlabeled_index)

view_num = size(X,1);
% for v=1:view_num
%     X{v} = X{v}(unlabeled_index,:);
% end
N = size(X{1},1);
e=7;
ker = struct('type','gauss','width',2); 

%初始化矩阵S
cell_s = cell(view_num,1); 
for t = 1:view_num
    cell_s{t,1} = zeros(N,N); 
end


%初始化矩阵D,L
cell_d = cell(view_num,1); 
cell_l = cell(view_num,1);
% tic;
%生成矩阵S及D
for t = 1:view_num
    acc_x = X{t,1};
    acc_s = cell_s{t,1};
    for i = 1:N
        x_vec = acc_x(i,:);
        x_dist = x_vec(ones(N,1),:);
        x_dist = sqrt( sum( (x_dist - acc_x).^2, 2 ) );
        [x_sort,index] = sort(x_dist);
        for e_nearest = 1:e
            j = index(e_nearest,:);
            acc_s(i,j) = kernel(ker,acc_x(i,:)',acc_x(j,:)');
        end
    end
    acc_d = diag( sum( acc_s));
    acc_l = zeros(size(acc_d));
    acc_l = acc_d - acc_s;
    cell_s{t,1} = acc_s;
    cell_d{t,1} = acc_d;
    cell_l{t,1} = acc_l;
end