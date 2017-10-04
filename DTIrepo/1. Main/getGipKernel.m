function krnl=getGipKernel(y)
% Adapted from code of:
% Twan van Laarhoven, Sander B. Nabuurs, Elena Marchiori,
% (2011) Gaussian interaction profile kernels for predicting drug–target interaction
% http://cs.ru.nl/~tvanlaarhoven/drugtarget2013/

	krnl = y*y';
	krnl = krnl / mean(diag(krnl));
	krnl = exp(-kernel_to_distance(krnl));

%     % an easier-to-read alternative
%     temp  = pdist(y,'euclidean');
%     temp  = temp .^ 2;
%     temp  = squareform(temp);
%     gamma = (sum(sum(y))) / size(y,1);
%     krnl  = exp(-temp./gamma);

end


function d=kernel_to_distance(k)
	di = diag(k);
	d = repmat(di,1,length(k)) + repmat(di',length(k),1) - 2 * k;
end