function [a] = readArray (filename)
fid = fopen (filename, 'r');
if (fid == -1) ; fprintf ('readArray: fopen failed!\n') ; return ; end
    
dmax = fscanf (fid, '%d');   % read number of dimensions
fgets (fid);                 % skip to next line
for d=1:dmax
    cmaxd(d) = fscanf (fid, '%d');  % read number of elements
    fgets (fid);                    % skip to next line
end
a = fscanf (fid, '%f');      % read in all elements
a = reshape (a, cmaxd);
fclose (fid);
endcon