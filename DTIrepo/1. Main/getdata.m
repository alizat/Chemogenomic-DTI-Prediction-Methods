function [Interaction,kCompound,kTarget,Did,Tid]=getdata(path,dataset)

    % ================================================================

    %------------------%
    % Adjacency Matrix %
    %------------------%

    % Extract the adjacency matrix...
    % The following command merges adjacency matrices that are 
    % extracted from multiple files which are specified using the 
    % pattern: [ path dataset '_admat_dgc.txt']
    newData1 = importdata ([ path dataset '_admat_dgc.txt']);
    Interaction = newData1.data;
    % Now, 'Interaction' has the matrix WITHOUT the column headers and
    % the row labels.

    % NOTE: The top-left-most cell in the extracted adjacency matrix 
    % is 'empty'

    % The uppermost row of the extracted adjacency matrix contains the 
    % column headers (drugs' IDs)
    Did = newData1.textdata(1,:);
    Did(1)=[];  % remove the first element (which is 'empty')

    % The leftmost row of the extracted adjacency matrix contains the 
    % row labels (targets' IDs)
    Tid=newData1.textdata(:,1);
    Tid(1)=[];  % remove the first element (which is 'empty')

    % ================================================================

    %----------------------------%
    % Compound Similarity Matrix %
    %----------------------------%

    % Extract the compound similarity matrices and merge them into one.
    newData1 = importdata ([ path dataset '_simmat_dc.txt']);
    kCompound = newData1.data;

    % ================================================================

    %--------------------------%
    % Target Similarity Matrix %
    %--------------------------%

    % Extract the target similarity matrices and merge them into one.
    newData1 = importdata ([ path dataset '_simmat_dg.txt']);
    kTarget = newData1.data;

    % ================================================================

    clear newData1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Checking for positive semi-definite. This is because, at present, these are "similarity" matrices rather
    %than "kernel" matrices. For now, we just add increments of a small epsilon to the diagonal, until the
    %matrix becomes positive semi-definite:

    %We had problems with kCompound having imaginary eigenvalues. Here, we symmetrize kCompound
    %before continuing:
    kCompound = (kCompound + kCompound') / 2;
    kTarget   = (kTarget   + kTarget')   / 2;

    epsilon = .1;
    while sum(eig(kCompound) >= 0) < size(kCompound,1) || isreal(eig(kCompound))==0
        kCompound = kCompound + epsilon*eye(size(kCompound,1));
    end
    while sum(eig(kTarget) >= 0) < size(kTarget,1) || isreal(eig(kTarget))==0
        kTarget = kTarget + epsilon*eye(size(kTarget,1));
    end