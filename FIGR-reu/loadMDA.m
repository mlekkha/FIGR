%
% Usage:  a = loadMDA ("myfile.mda")
% Loads a multidimensional array from a MDA-format file.
% Example of a MDA file:
%
% 3 dims
% 4 elems
% 2 elems
% 3 elems
% 0. 1. 2. 3.
% 4. 5. 6. 7.0	
%
% 0. 1. 2. 3.
% 4. 5. 6. 7.0	
%
% 0. 1. 2. 3.
% 4. 5. 6. 7.0	

function [a] = loadMDA (filename)
fid = fopen (filename, 'r');
if (fid == -1) ; fprintf (2,'loadMDA: fopen failed!\n') ; return ; end
        
dmax = fscanf (fid, '%d');   % read number of dimensions
fgets (fid);                 % skip to next line
for d=1:dmax
    cmaxd(d) = fscanf (fid, '%d');  % read number of elements
    fgets (fid);                    % skip to next line
end
a = fscanf (fid, '%f');      % read in all elements
a = reshape (a, cmaxd);
fclose (fid);
fprintf ('loadMDA: read array of size [%s] from file %s\n', num2str(cmaxd), filename);
end