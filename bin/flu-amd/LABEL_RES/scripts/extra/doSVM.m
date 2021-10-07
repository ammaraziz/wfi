#!/opt/local/bin/octave -qf
% Samuel Shepard - doSVM - 7.2011
% Tests unknown sequences for the genotyping system.
% The load_matrix function is given below and taken
% from the SHOGUN toolkit help files.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Courtesy of the SHOGUN Toolkit.
function [matrix] = load_matrix(fname)
%function [matrix] = load_matrix(fname)
%
% load a matrix from file
%
% this implementation is prone to become rather slow with large data,
% but our expected toy data here is not large. :)
	fid=fopen(fname, 'r');
	matrix=[];
	while 1
		line=fgetl(fid);
		if ~ischar(line),   break,   end

		converted=str2num(line);
		if isempty(converted)
			matrix=[matrix, line'];
		else
			matrix=[matrix; converted];
		end
	end
	fclose(fid);
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


size_cache=10;
use_bias=false;
epsilon=1e-5;
width=2.1;
degree=20;

%poly
inhomogene=true;
use_normalization=true;

alist=argv();
path=sprintf('%s/resources/training_data/%s_', alist{3}, alist{1});
labels=load_matrix( [path, 'training_labels.dat'] );
traindata=load_matrix( [path,'training.dat'] );
input=sprintf( '%s/resources/test_data/%s/%s_%s.dat', alist{3}, alist{2}, alist{1}, alist{2} );
testdata=load_matrix( input );

sg('new_classifier', 'GMNPSVM');
sg('set_kernel', 'POLY', 'REAL', size_cache, degree, inhomogene, use_normalization);
sg('set_features', 'TRAIN', traindata);
sg('set_labels', 'TRAIN', labels);
sg('train_classifier');

sg('set_features', 'TEST', testdata);
R=sg('classify');

L = length(R);

fid = fopen( sprintf('%s/resources/test_data/%s/%s_predictions.tab', alist{3}, alist{2},alist{2}), "a");
for i=1:1:L
	fprintf(fid,"%d\t%d\t%s\n",i,R(i),alist{1});
end
fclose(fid);
