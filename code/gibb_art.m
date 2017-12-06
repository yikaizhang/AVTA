function [dataset]=gibb_art(A,K,word_num_m,sample_num)
%%********Generate semi synthetic data******
%Input:
%   A: word topic distribution
%   word num_m: Mean length of the documents
%   sample_num: The number of samples to be generated
%Output:n
%   dataset: The generated dataset. A column of 'dataset' is a document.

ALPHA=0.01;

word_num=word_num_m;

dirichlet_val=dirichletrnd(ALPHA*ones(K,1), sample_num);

dataset=mnrnd(word_num,(A*dirichlet_val)');
dataset=dataset';