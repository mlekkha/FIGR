% At the MATLAB command line, just cd to the folder containing this script,
% and type
%
%   readGeneExprFiles
%
% Conventions:
%
%   t    timepoint index (3) or time value (3.75 minutes)
%   n    nucleus index (2) or nucleus label (8227)
%   s    sample index (3) or sample ID ("Sample_4_8226")
%   u    unique-datapoint index (ideally umax=nmax*tmax)
%   g    gene index


% EVENTUALLY WE WANT TO RETURN THESE:
% [tt, nn, x_ntg, x_ntgMean] = 
function [xntg tt geneNames] = readGeneExprFiles ()

fnameTN = "tn.txt";
fnameXNTG = "xntg.txt";
%fnameTN = "simple_tn.txt";
%fnameXNTG = "simple_xntg.txt";

%======== READ THREE-COLUMN FILE CONTAINING SAMPLE/TIME/NUCLEUS INFO
fid = fopen (fnameTN, 'r');        %Opens File containing time points
matrix = textscan (fid, '%s%s%s', 'Delimiter', '\t', 'HeaderLines', 1, 'MultipleDelimsAsOne', 1);
fclose (fid);
ss = strtrim (matrix{1,1});      % list of sample IDs for each sample s (whitespace removed!)
ts = str2double (matrix{1,2});   % list of times t_s for each sample s
ns = str2double (matrix{1,3});   % list of nuclei n_s for each sample s
tt = unique (ts);  % actual time T_t corresponding to timepoint index t
nn = unique (ns);  % nucleus label N_n corresponding to nucleus index n
smax = numel (ss); % number of samples
tmax = numel (tt); % number of unique timepoints
nmax = numel (nn); % number of unique nuclei

%======== COMPARE SAMPLE NAMES BETWEEN TWO FILES
fid = fopen(fnameXNTG , 'r');
dummy = fgetl (fid); 
dummy = split (strtrim (dummy)); % split into strings
ssFromMainFile = dummy(2:end);
assert (numel(ssFromMainFile) == smax, 'ERROR: incompatible numbers of samples in data files!');
assert (isequal (ssFromMainFile, ss), 'ERROR: differing sample names in data files!');
%======== READ EACH LINE OF THE FILE AND PARSE IT INTO THE ARRAY
geneNames = [];                          % could be "gg"
x_gs = [];
while true
    dummy = fgetl (fid);                 % read one line from the file (as char vec)
    if ~ischar(dummy); break; end
    dummy = split (strtrim (dummy));     % split line into strings
    geneNames = [geneNames ; dummy(1)];  % first item is gene name
    dummy = str2double ( dummy(2:end) )'; % remaining items are gene expression levels
    x_gs = [x_gs ; dummy];
end
fclose (fid);
gmax = numel (geneNames);

%======== AVERAGE OVER REPLICATES
[tnu, ia, ic] = unique ([ts ns], 'rows', 'stable');
umax = numel (ia);  % number of unique datapoints
x_gu = zeros ([gmax umax]);
for g = 1:gmax
	x_gu(g,:) = accumarray (ic, x_gs(g,:), [], @mean)';
end

%======== RESHAPE x(n,t,g) INTO x(g,u)
xntg = zeros ([nmax tmax gmax]);
for u = 1:umax
    tn = tnu(u,:);  % time and position corresponding to sample s
    t = tn(1);
    n = tn(2);
    indext = find (tt==t);  % time index corresponding to sample s
    indexn = find (nn==n);
    xntg (indexn,indext,:) = x_gu (:,u);
end	

%======== DUMP DEBUGGING INFO
fprintf ("\n", smax);
fprintf ("Number of samples     smax = %d\n", smax);
fprintf ("Number of uniques     umax = %d\n", umax);
fprintf ("Number of nuclei      nmax = %d\n", nmax);
fprintf ("Number of timepoints  tmax = %d\n", tmax);
fprintf ("Number of genes       gmax = %d\n", gmax);

format bank;
disp ("Before averaging over replicates:");
%disp ("x_gs = "); disp (x_gs);
%fprintf ("Replica indices: "); fprintf ("%d ", ic); fprintf ("\n");
%fprintf ("Replica counts:  "); fprintf ("%d ", replicates); fprintf ("\n");
disp ("After averaging over replicates:");
%disp ("x_gu = "); disp (x_gu);
disp ("Parsed into FIGR format:");
%disp ("x_ntg = "); disp (xntg);
format;
end





%mapReplicates = replicates (ic);
%result = [ts ns mapReplicates];  % 

%dummy = textscan (fid, fmtstr, 'MultipleDelimsAsOne', 1);
% textscan (fid, '%s', 1, 'Delimiter', '\n'); 
%dummy = dummy{1};
