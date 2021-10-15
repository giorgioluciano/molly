clc; clear; close all;


 EdgeTable = table([1 2; 2 3; 3 5; 4 5; 6 7],'VariableNames',{'EndNodes'});
 G = graph(EdgeTable);
 plot(G)
 arcs = [1 2; 2 3; 3 5; 4 5; 6 7];
 cc=conncomp(graph(table(arcs,'VariableNames',{'EndNodes'})));


