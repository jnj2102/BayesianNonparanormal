%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to calculate the mean and standard errors of the simulation
%combinations for the Bayesian nonparanormal method
% Author: Jami Jackson Mulgrave
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%tenpercent model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%p=25 n=50

load('BayesNonpar_Simulation_n50_p25_tenpercent_SpikeSlab_1to100_final.mat');

%Vertically concatenate the row vectors of the Bayes L1 loss 
SE_n50_p25_tenpercent_vector = zeros([reps,1]); 
SP_n50_p25_tenpercent_vector = zeros([reps,1]); 
MCC_n50_p25_tenpercent_vector = zeros([reps,1]);
L1_n50_p25_tenpercent_vector = zeros([reps,1]);
total_time_n50_p25_tenpercent_vector = zeros([reps,1]);



for iters = 1:reps
    total_time_n50_p25_tenpercent_vector(iters) = total_time_n50_p25_tenpercent_finalanalysis(iters);
    L1_n50_p25_tenpercent_vector(iters) = L1_finalanalysis(iters);
    SP_n50_p25_tenpercent_vector(iters) = SP_matrix_finalanalysis(iters);
	SE_n50_p25_tenpercent_vector(iters) =  SE_matrix_finalanalysis(iters);
	MCC_n50_p25_tenpercent_vector(iters) =   MCC_matrix_finalanalysis(iters);

end

%findd the mean time and L1

total_time_n50_p25_tenpercent_mean_rounded = round2(mean(total_time_n50_p25_tenpercent_vector),.01);
total_time_n50_p25_tenpercent_SE_rounded = round2(std(total_time_n50_p25_tenpercent_vector)/sqrt(length(total_time_n50_p25_tenpercent_vector)),.01);


L1_n50_p25_tenpercent_mean_rounded = round2(mean(L1_n50_p25_tenpercent_vector),.01);
L1_n50_p25_tenpercent_SE_rounded = round2(std(L1_n50_p25_tenpercent_vector)/sqrt(length(L1_n50_p25_tenpercent_vector)),.01);



%find the mean and SE for  SE, SP, MCC

%Since the cluster has an older version of Matlab, the round2 function is different.  I have to use round2 and use .01 for the number of significant digits

SE_n50_p25_tenpercent_mean_rounded = round2(mean(SE_n50_p25_tenpercent_vector),.01);
SE_n50_p25_tenpercent_SE_rounded = round2(std(SE_n50_p25_tenpercent_vector)/sqrt(length(SE_n50_p25_tenpercent_vector)),.01);

SP_n50_p25_tenpercent_mean_rounded = round2(mean(SP_n50_p25_tenpercent_vector),.01);
SP_n50_p25_tenpercent_SE_rounded = round2(std(SP_n50_p25_tenpercent_vector)/sqrt(length(SP_n50_p25_tenpercent_vector)),.01);

%have to exclude NaNs when reporting MCC.
MCC_n50_p25_tenpercent_mean_rounded = round2(nanmean(MCC_n50_p25_tenpercent_vector),.01);
MCC_n50_p25_tenpercent_SE_rounded = round2(nanstd(MCC_n50_p25_tenpercent_vector)/sqrt(length(MCC_n50_p25_tenpercent_vector(~isnan(MCC_n50_p25_tenpercent_vector)))),.01);


save('Bsplines_paper_Calcs_p25_n50_tenpercent_SpikeSlab_1to100.mat',...
    'total_time_n50_p25_tenpercent_mean_rounded',... 
'total_time_n50_p25_tenpercent_SE_rounded',...
'L1_n50_p25_tenpercent_mean_rounded',... 
'L1_n50_p25_tenpercent_SE_rounded',...
   'SE_n50_p25_tenpercent_mean_rounded',...
'SE_n50_p25_tenpercent_SE_rounded', 'SP_n50_p25_tenpercent_mean_rounded', 'SP_n50_p25_tenpercent_SE_rounded', 'MCC_n50_p25_tenpercent_mean_rounded',...
 'MCC_n50_p25_tenpercent_SE_rounded');

