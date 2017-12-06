function nodes_from_mat_squares(~)
%IMPORTFILE(FILETOREAD1)
%  Imports data from the specified file
%
%name='15x15';
% First going to make the mesh_node file
x = load('/Users/crhea/Documents/Grad/Math/530_FEM/Final/Code/mesh/matlab/fem_mini.mat','x0');
y = load('/Users/crhea/Documents/Grad/Math/530_FEM/Final/Code/mesh/matlab/fem_mini.mat','y0');
x = x.('x0');
y = y.('y0');
% Create new variables in the base workspace from those fields.
FileID = fopen('/Users/crhea/Documents/Grad/Math/530_FEM/Final/Code/mesh/mesh_nodes_mini.txt','w');
for i = 1:length(x)
    fprintf(FileID,'%6f ',x(i));
    fprintf(FileID,'%6f\n',y(i));
end

%now for the connectivity file
blk = load('/Users/crhea/Documents/Grad/Math/530_FEM/Final/Code/mesh/matlab/fem_mini.mat','blk01');
blk = blk.('blk01');
fileid = fopen('/Users/crhea/Documents/Grad/Math/530_FEM/Final/Code/mesh/mesh_connectivity_mini.txt','w');
disp(length(blk))
for i = 1:length(blk)
    fprintf(fileid,'%d ',blk(i*4-3));
    fprintf(fileid,'%d ',blk(i*4-2));
    fprintf(fileid,'%d ',blk(i*4-1));
    fprintf(fileid,'%d \n',blk(i*4));
end

%now lets get our boundaries
top = load('/Users/crhea/Documents/Grad/Math/530_FEM/Final/Code/mesh/matlab/fem_mini.mat','sselem01');
top = top.('sselem01');
left = load('/Users/crhea/Documents/Grad/Math/530_FEM/Final/Code/mesh/matlab/fem_mini.mat','sselem02');
left = left.('sselem02');
bottom = load('/Users/crhea/Documents/Grad/Math/530_FEM/Final/Code/mesh/matlab/fem_mini.mat','sselem03');
bottom = bottom.('sselem03');
right = load('/Users/crhea/Documents/Grad/Math/530_FEM/Final/Code/mesh/matlab/fem_mini.mat','sselem04');
right = right.('sselem04');

fileid2 = fopen('/Users/crhea/Documents/Grad/Math/530_FEM/Final/Code/mesh/sides_mini.txt','w');

for i = 1:length(top)
    fprintf(fileid2,'%d ',top(i));
    fprintf(fileid2,'%d ',left(i));
    fprintf(fileid2,'%d ',bottom(i));
    fprintf(fileid2,'%d \n',right(i));
end
fprintf(fileid2,'\n');



end


