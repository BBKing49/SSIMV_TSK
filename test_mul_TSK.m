function [Y_te] = test_mul_TSK( data ,TSK_cell,  view_nums, c, unlabeled_index)
%view_nums视角数，c为类别数

for t = 1:view_nums
    test_cell{t,1} = data{t,1}(unlabeled_index,:);
end
Y_te = zeros( size( test_cell{1,1},1), c);

for t = 1:view_nums
    acc_pg = TSK_cell{t,1};
    acc_v = TSK_cell{t,2};
    acc_b = TSK_cell{t,3};
    acc_w = TSK_cell{t,4};
    acc_x = test_cell{t,1};
%     acc_x = fromXtoZ(acc_x,acc_v,acc_b);
    FS_output = acc_x*acc_pg;
    Y_te = Y_te + acc_w*FS_output;
end