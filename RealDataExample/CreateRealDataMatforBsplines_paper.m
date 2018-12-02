%Code to convert the real data to a mat file for further use

data = readtable('gb-2004-5-11-r92-s1_removedrows.txt'); % read in the data

data_table = data(:, 7:124);

data_matrix = table2array(data_table);


save('Bsplines_paper_realdata_initialdata.mat', '-v7.3');