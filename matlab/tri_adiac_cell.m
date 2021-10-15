function [ad_tri,ext_mat_tri]=tri_adiac_cell(mat_ver,mat_tri)

num_ver=size(mat_ver,1); num_tri=size(mat_tri,1);

succ=[2 3 1]; t=1:num_tri;

t_t=spalloc(num_ver,num_ver,3*num_tri);
for i=1:3
    t_t=t_t+sparse(mat_tri(:,i),mat_tri(:,succ(i)),t,num_ver,num_ver,num_tri);
end

for i=1:3
    index_edge=sub2ind(size(t_t),mat_tri(:,succ(i)),mat_tri(:,i));
    mat_tri(:,i+3)=t_t(index_edge);
end
ext_mat_tri=mat_tri;
mat_tri=mat_tri(:,4:6);
%ad_tri=cell(1,size(mat_tri,1)); %if you want a cell array as output
ad_tri=zeros(size(mat_tri));

for i=1:size(mat_tri,1)
    % FIX n1
    bb=[];
    cc=0;
    for ii=1:3
        if mat_tri(i,ii)~=0
            cc=cc+1;
            bb(1,cc)=mat_tri(i,ii);
        end
    end
    cc=0;
    %ad_tri{1,i}=bb; %if you want a cell array as output
    ad_tri(i,:)=bb;
end