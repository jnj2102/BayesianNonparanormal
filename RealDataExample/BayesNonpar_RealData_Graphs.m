%Code to create the estimated graphs for the Real Data application
%
%Author: Jami Jackson Mulgrave
clear;

addpath('C:/Users/jnjac/Documents/MATLAB/paul-kassebaum-mathworks-circularGraph-3a7926b');

load('RealData_Frequentist.mat');
load('BayesNonpar_RealData.mat');
load('RealData_BDGraph.mat');

iters = 1;
edge_matrix_finalanalysis = edge_matrix_finalanalysis{iters};

%Bayesian nonparanormal graph
%edge_matrix_finalanalysis(logical(eye(size(edge_matrix_finalanalysis)))) = 0;

edge_matrix_finalanalysis_matrix = double(edge_matrix_finalanalysis);

%call the circular graph

figure;

myColorMap = repmat([0 0 0], [p,1]);

H = circularGraph(edge_matrix_finalanalysis_matrix, 'Colormap',myColorMap);


savefig('RealData_circularGraph_Bayes.fig')

figs = openfig('RealData_circularGraph_Bayes.fig');
   saveas(figs, 'RealData_circularGraph_Bayes.png'); %MATLAB doesn't recognize 
   %circular graph in the saveas function, so I need to save it as a fig
   %and then convert it to another extension that can be opened elsewhere



%Read in the frequentist STARS nonpararnormal model

figure;

stars_graph = circularGraph(edge_matrix_est_stars, 'Colormap',myColorMap);

savefig('RealData_circularGraph_Frequentist_stars.fig')

figs = openfig('RealData_circularGraph_Frequentist_stars.fig');
   saveas(figs, 'RealData_circularGraph_Frequentist_stars.png');
   

%Read in the frequentist RIC nonpararnormal model

figure;

ric_graph = circularGraph(edge_matrix_est_ric, 'Colormap',myColorMap);

savefig('RealData_circularGraph_Frequentist_ric.fig')

figs = openfig('RealData_circularGraph_Frequentist_ric.fig');
   saveas(figs, 'RealData_circularGraph_Frequentist_ric.png');

%Read in the frequentist EBIC nonpararnormal model


figure;

ebic_graph = circularGraph(edge_matrix_est_ebic, 'Colormap',myColorMap);

savefig('RealData_circularGraph_Frequentist_ebic.fig')

figs = openfig('RealData_circularGraph_Frequentist_ebic.fig');
   saveas(figs, 'RealData_circularGraph_Frequentist_ebic.png');%How many edges do we get for each?

   
   %BDGraph
   
   
   figure;

BDGraph_graph = circularGraph(edge_matrix_selected, 'Colormap',myColorMap);

savefig('RealData_circularGraph_BDGraph.fig')

figs = openfig('RealData_circularGraph_BDGraph.fig');
   saveas(figs, 'RealData_circularGraph_BDGraph.png');
   
   %How many different edges per graph?
   

indmx = reshape(1:p^2,p,p); 
  upperind = indmx(triu(indmx,1)>0);  %do not include the diagonal
 
  sum_edges_bayes =  sum(edge_matrix_finalanalysis(upperind) == 1);

    sum_edges_stars =  sum(edge_matrix_est_stars(upperind) == 1);

        sum_edges_ric =  sum(edge_matrix_est_ric(upperind) == 1);

          sum_edges_ebic =  sum(edge_matrix_est_ebic(upperind) == 1);

                                        sum_edges_BDGraph =  sum(edge_matrix_selected(upperind) == 1);
