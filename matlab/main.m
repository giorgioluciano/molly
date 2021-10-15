clc; clear; close all;

T = load("sprout_tetrahedra_atoms_32.txt");

pre_process = 0;
plot_tetra = 0;
K = 25;


%% PRE-PROCESSING

if pre_process==1
    
    % 1)Read the number of vertices and facets:
    fid = fopen("triangulatedSurf.txt");
    linenum = 4;    %The file MUST have the number of vertices in the 4th row
    C = textscan(fid,'%f',3,'delimiter','\n', 'headerlines',linenum-1);
    if (floor(C{1}(2))==C{1}(2)) %meaning, we have the number of vertices!
        num_vert = C{1}(1);
        num_facets = C{1}(2);
    else
        error("UNEXPECTED FLOAT...");
    end
    fclose(fid);
    
    % 2) Extract vertices and facets:
    system(strcat("sed -n '", num2str(5), ",", num2str(4+num_vert), " p'", " ./", "triangulatedSurf.txt", " > ", "./vertices.txt"));
    system(strcat("sed -n '", num2str(5+num_vert), ",", num2str(4+num_vert+num_facets), " p'", " ./", "triangulatedSurf.txt", " > ", "./facets.txt"));
    
end

% Load vertices and facets:
V = load("vertices.txt");
F = load("facets.txt");


%% BARYCENTER OF TETRAHEDRA (OF THE CCs) AND OF FACETS (OF THE SES)
B_tet_CC = (T(:,1:3)+T(:,4:6)+T(:,7:9)+T(:,10:12))/4;
B_fac_SES = (V(F(:,2)+1,:)+V(F(:,3)+1,:)+V(F(:,4)+1,:))/3;

%Point-triangle distance:
%https://it.mathworks.com/matlabcentral/fileexchange/22857-distance-between-a-point-and-a-triangle-in-3d
[Idx_SES2TET, D_SES2TET] = knnsearch(B_tet_CC, B_fac_SES,'K',K);



%% FACETS OF THE CCs
%https://doc.cgal.org/latest/TDS_3/index.html

% First facet:
F1 = T(:, [1,2,3,4,5,6,10,11,12]);
% Second facet:
F2 = T(:, [4,5,6,7,8,9,10,11,12]);
% Third facet:
F3 = T(:, [1,2,3,7,8,9,4,5,6]);
% Fourth facet:
F4 = T(:, [1,2,3,10,11,12,7,8,9]);

% Normals to facets:
N_facets_F1 = cross(F1(:,4:6)-F1(:,1:3), F1(:,7:9)-F1(:,1:3));
N_facets_F2 = cross(F2(:,4:6)-F2(:,1:3), F2(:,7:9)-F2(:,1:3));
N_facets_F3 = cross(F3(:,4:6)-F3(:,1:3), F3(:,7:9)-F3(:,1:3));
N_facets_F4 = cross(F4(:,4:6)-F4(:,1:3), F4(:,7:9)-F4(:,1:3));

n = vecnorm(N_facets_F1,2,2); N_facets_F1 = N_facets_F1./n;
n = vecnorm(N_facets_F2,2,2); N_facets_F2 = N_facets_F2./n;
n = vecnorm(N_facets_F3,2,2); N_facets_F3 = N_facets_F3./n;
n = vecnorm(N_facets_F4,2,2); N_facets_F4 = N_facets_F4./n;

clear n;

%Baricenters of the facets:
B_F1 = [(F1(:,1)+F1(:,4)+F1(:,7))/3,(F1(:,2)+F1(:,5)+F1(:,8))/3,(F1(:,3)+F1(:,6)+F1(:,9))/3]; B_F1 = round(B_F1, 3);
B_F2 = [(F2(:,1)+F2(:,4)+F2(:,7))/3,(F2(:,2)+F2(:,5)+F2(:,8))/3,(F2(:,3)+F2(:,6)+F2(:,9))/3]; B_F2 = round(B_F2, 3);
B_F3 = [(F3(:,1)+F3(:,4)+F3(:,7))/3,(F3(:,2)+F3(:,5)+F3(:,8))/3,(F3(:,3)+F3(:,6)+F3(:,9))/3]; B_F3 = round(B_F3, 3);
B_F4 = [(F4(:,1)+F4(:,4)+F4(:,7))/3,(F4(:,2)+F4(:,5)+F4(:,8))/3,(F4(:,3)+F4(:,6)+F4(:,9))/3]; B_F4 = round(B_F4, 3);

if plot_tetra==1
    for i=1:size(B_tet_CC,1)
        hold on;
        DT = delaunayTriangulation([T(i,1); T(i,4); T(i,7); T(i,10)],...
            [T(i,2); T(i,5); T(i,8); T(i,11)], [T(i,3); T(i,6); T(i,9); T(i,12)]);
        tetramesh(DT, 'FaceAlpha', 0.2);
        %textscatter3([T(i,1); T(i,4); T(i,7); T(i,10)],...
        %    [T(i,2); T(i,5); T(i,8); T(i,11)], [T(i,3); T(i,6); T(i,9); T(i,12)], ...
        %    string(0:3));
        
        quiver3(B_F1(i,1), B_F1(i,2), B_F1(i,3), N_facets_F1(i,1), N_facets_F1(i,2), N_facets_F1(i,3));
        quiver3(B_F2(i,1), B_F2(i,2), B_F2(i,3), N_facets_F2(i,1), N_facets_F2(i,2), N_facets_F2(i,3));
        quiver3(B_F3(i,1), B_F3(i,2), B_F3(i,3), N_facets_F3(i,1), N_facets_F3(i,2), N_facets_F3(i,3));
        quiver3(B_F4(i,1), B_F4(i,2), B_F4(i,3), N_facets_F4(i,1), N_facets_F4(i,2), N_facets_F4(i,3));
        axis equal;
        disp("");
    end
end



%% COLORING POCKETS

%Coloring what is inside the current CC:
tic
colors = zeros(size(B_fac_SES,1),1);
for k=1:size(B_fac_SES,1)
    for l=1:K
        if ...
                (B_fac_SES(k,:)-B_F1(Idx_SES2TET(k,l),:))*N_facets_F1(Idx_SES2TET(k,l),:)'<0 && ...
                (B_fac_SES(k,:)-B_F2(Idx_SES2TET(k,l),:))*N_facets_F2(Idx_SES2TET(k,l),:)'<0 && ...
                (B_fac_SES(k,:)-B_F3(Idx_SES2TET(k,l),:))*N_facets_F3(Idx_SES2TET(k,l),:)'<0 && ...
                (B_fac_SES(k,:)-B_F4(Idx_SES2TET(k,l),:))*N_facets_F4(Idx_SES2TET(k,l),:)'<0
            colors(k)=1;
            break
        end
    end
end
toc


%Filling holes:

edges = [];

tic
[ad_tri,ext_mat_tri]=tri_adiac_cell(V,F(:,2:4)+1);
col_ad_tri = colors(ad_tri);

% First facet
a1= find(and(any(colors == 0, 2), any(col_ad_tri(:,1) == 0, 2)));
a2 = find(and(any(colors == 1, 2), any(col_ad_tri(:,1) == 1, 2)));
edges = [edges; a1 ad_tri(a1,1)];
edges = [edges; a2 ad_tri(a2,1)];

% Second facet:
a1= find(and(any(colors == 0, 2), any(col_ad_tri(:,2) == 0, 2)));
a2 = find(and(any(colors == 1, 2), any(col_ad_tri(:,2) == 1, 2)));
edges = [edges; a1 ad_tri(a1,2)];
edges = [edges; a2 ad_tri(a2,2)];

% Third facet:
a1= find(and(any(colors == 0, 2), any(col_ad_tri(:,3) == 0, 2)));
a2 = find(and(any(colors == 1, 2), any(col_ad_tri(:,3) == 1, 2)));
edges = [edges; a1 ad_tri(a1,3)];
edges = [edges; a2 ad_tri(a2,3)];
toc

cc=conncomp(graph(table(edges,'VariableNames',{'EndNodes'})));
count_cc = zeros(max(cc),1);
for i=1:max(cc)
    count_cc(i)=sum(cc==i);
end
id_max_count_cc = find(count_cc==max(count_cc));

for i=1:max(cc)
    if i~=id_max_count_cc
        colors(cc==i)=1;
    end
end

facets_gray = colors==0;
facets_red = colors==1;

F_colored = [F, ones(size(F,1),4)];
F_colored(facets_gray,5:7) = 200*F_colored(facets_gray,5:7);
F_colored(facets_red,5) = 255*F_colored(facets_red,5);
F_colored(facets_red,6) = 0*F_colored(facets_red,6);
F_colored(facets_red,7) = 0*F_colored(facets_red,7);
F_colored(:,end) = 255*F_colored(:,end);


%% GETTING RID OF UNINTERESTING HOLES



%Saving COFF:
fileID = fopen('colored_triangulation_32.txt','w');
fprintf(fileID,'COFF\n');
fprintf(fileID, '%d %d 0\n', size(V,1), size(F,1));
for i=1:size(V,1)
    fprintf(fileID, '%f %f %f\n', V(i,1), V(i,2), V(i,3));
end
for i=1:size(F_colored,1)
    for j=1:size(F_colored,2)-1
        fprintf(fileID, '%d ', F_colored(i,j));
    end
    fprintf(fileID, '%d\n', F_colored(i,size(F_colored,2)));
end

