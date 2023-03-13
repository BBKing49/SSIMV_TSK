function [Y] = test_TSK_FS( test_data , pg, v, b)
%test_dataÎªN*D£¬Êä³öY£ºN*C 
test_data = fromXtoZ(test_data,v,b);
Y = test_data*pg;
