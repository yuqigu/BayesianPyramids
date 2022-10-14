% This script processes the splice junction data before using the CORELS
% package in Python
Y_with_type = csvread('slice_no_n.csv');
Y_mat_d4 = Y_with_type(:, 2:end);
seq_type = Y_with_type(:, 1);

splice_Y1 = [Y_mat_d4, seq_type==1];
splice_Y2 = [Y_mat_d4, seq_type==2];
splice_Y3 = [Y_mat_d4, seq_type==3];

csvwrite('splice_Y1.csv', splice_Y1);
csvwrite('splice_Y2.csv', splice_Y2);
csvwrite('splice_Y3.csv', splice_Y3);


% -- binary data -- %
Y_bin_manycol = [Y_mat_d4==1, Y_mat_d4==2, Y_mat_d4==3];
splice_bin_Y1 = [Y_bin_manycol, seq_type==1];
splice_bin_Y2 = [Y_bin_manycol, seq_type==2];
splice_bin_Y3 = [Y_bin_manycol, seq_type==3];

csvwrite('splice_bin_Y1.csv', splice_bin_Y1);
csvwrite('splice_bin_Y2.csv', splice_bin_Y2);
csvwrite('splice_bin_Y3.csv', splice_bin_Y3);