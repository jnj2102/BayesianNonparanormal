%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boxplots for the BayesNonpar paper
% Author: Jami Jackson Mulgrave
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%twopercent p100 n500 SpikeSlab 

load('BayesNonpar_p100_n500_twopercent_SpikeSlab_final.mat');


SP_value = SP_matrix_finalanalysis';
SP_Type = repmat({'Specificity'}, [reps,1]);
SE_value = SE_matrix_finalanalysis';
SE_Type = repmat({'Sensitivity'}, [reps,1]);
MCC_value = MCC_matrix_finalanalysis';
MCC_Type = repmat({'MCC'}, [reps,1]);

Edges_value = [SP_value; SE_value; MCC_value];
Edges_Type = [SP_Type; SE_Type; MCC_Type];

[number_elements, ~] = size(Edges_value);

Sparsity = repmat({'Percent'}, [number_elements,1]);
Method = repmat({'Spike Slab'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n500_twopercent_SpikeSlab = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [table_p100_n500_twopercent_SpikeSlab];


clearvars -except combine_tables


%circle p100 n500 SpikeSlab 

load('BayesNonpar_p100_n500_circle_SpikeSlab_final.mat');

SP_value = SP_matrix_finalanalysis';
SP_Type = repmat({'Specificity'}, [reps,1]);
SE_value = SE_matrix_finalanalysis';
SE_Type = repmat({'Sensitivity'}, [reps,1]);
MCC_value = MCC_matrix_finalanalysis';
MCC_Type = repmat({'MCC'}, [reps,1]);

Edges_value = [SP_value; SE_value; MCC_value];
Edges_Type = [SP_Type; SE_Type; MCC_Type];

[number_elements, ~] = size(Edges_value);

Sparsity = repmat({'Circle'}, [number_elements,1]);
Method = repmat({'Spike Slab'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n500_circle_SpikeSlab = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n500_circle_SpikeSlab];


clearvars -except combine_tables

%AR1 p100 n500 SpikeSlab 

load('BayesNonpar_p100_n500_AR1_SpikeSlab_final.mat');

SP_value = SP_matrix_finalanalysis';
SP_Type = repmat({'Specificity'}, [reps,1]);
SE_value = SE_matrix_finalanalysis';
SE_Type = repmat({'Sensitivity'}, [reps,1]);
MCC_value = MCC_matrix_finalanalysis';
MCC_Type = repmat({'MCC'}, [reps,1]);

Edges_value = [SP_value; SE_value; MCC_value];
Edges_Type = [SP_Type; SE_Type; MCC_Type];

[number_elements, ~] = size(Edges_value);

Sparsity = repmat({'AR1'}, [number_elements,1]);
Method = repmat({'Spike Slab'}, [number_elements,1]);
Dimension = repmat(p, [number_elements,1]);

table_p100_n500_AR1_SpikeSlab = table(Edges_value, Method, Edges_Type, Sparsity, Dimension);

combine_tables = [combine_tables; table_p100_n500_AR1_SpikeSlab];



%Try writing table as a csv file to read into R

writetable(combine_tables,'BayesNonpar_Boxplot_Edges_p100_n500.csv') 

