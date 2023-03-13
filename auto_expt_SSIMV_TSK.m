clear;
clc;

data_files = {'Dermatology'};

tic;
for data_file = data_files
    load(['data/' data_file{1} '.mat']);
     labeled_rate = 0.3;
    for incomplete_rate = [ 0.5]
        for i = 1:3
            output_file = ['results/' data_file{1} int2str(i) '_' num2str(incomplete_rate)  '.mat'];  
            [best_result,TSK_result ] = expt_SSIMV_TSK( data, labels, incomplete_rate, labeled_rate, output_file );
            save(output_file, 'best_result');
        end
    end
end
toc;