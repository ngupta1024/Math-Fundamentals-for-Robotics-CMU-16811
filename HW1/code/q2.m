clear;

% Test cases
A=[[1,-1,0];[0,2,-1]; [1,0,-1/2]];
%A=[[-1,1,0,0];[-1,0,1,0];[0,-4,1,0];[0,-1,0,1];[0,0,-2,1]];
%A=[[2,2,5];[1,1,5];[3,2,5]];
%A=[[-3,-4,-1];[2,3,1];[3,5,2]];
%A=[[1,-1,0];[0,1,-2];[1,0,-2]];

% Decompose A with PA=LDU
[P,L,D,U]=LDU_Decompose(A);
% Decompose A with A=USV
[U2,S,V]=SVD_Decompose(A);

% Check the correctness of the LUD deposition
status1=check_LDU(P,A,L,D,U);
% Check the correctness of the SVD deposition
status2=check_SVD(A,U2,S,V);


function [P,L,D,U]=LDU_Decompose(A)
% Inputs: 
%   A: The input m*n matrix that needs LDU decomposition
%
% Outputs:
%   P: The permutation matrix that is used for shifting the rows of A
%   L: The decomposition result L, which is a lower triangular matrix
%   D: The decomposition result D, which is a diagonal matrix
%   U: The decomposition result U, which is a upper triangular matrix
%      The returned matrices make PA = LDU
    [m,n]=size(A);
    A_=A;
    P=eye(m);
    L=eye(m);
    for k=1:min(m,n)
        if A_(k,k)==0
            % need to swap the rows to make A_(k,k) is non-zero
            row_to_cahnge=k;
            for t=k+1:m
                if A_(t,k)~=0
                    row_to_cahnge=t;
                    break
                end
            end
            [P,A_,L]=row_change(A_,P,L,row_to_cahnge,k,true);
        end
        % process all the other rows with gaussian elimination method
        for t=k+1:m
            [P,A_,L]=row_change(A_,P,L,t,k,false);
        end
    end
    % calculate D and U since DU=A_
    D=zeros(m,m);
    diagA_=diag(diag(A_(1:min(m,n),1:min(m,n))));
    r=rank(diagA_);
    D(1:r,1:r)=diagA_(1:r,1:r);
    U=eye(m,n);
    U(1:r,:)=D(1:r,1:r)\A_(1:r,:);
end

function [P_,A_,L_]=row_change(A_,P_,L_,row2,row1,isSwap)
% Inputs:
%   A_: The square and invertible matrix that needs to shift the rows
%   P_: The permutation matrix that used for shifting the rows of A_
%   L_: The matrix that makes L_*A_=A, it is a lower triangular matrix
%   row2: The index of target row to be processed with guassian elimination
%   row1: The index of base row used for processing row2
%   isSwap: Be in type of bool. If true, then just swap row1 and row2
% Outputs:
%   P_: The permutation matrix after processing row2
%   A_: The matrix after applying guassian elimination to row2
%   L_: The corresponding matrix that makes L_*A_=A after row processing
    if isSwap
        % operations when need to swap two rows
        A_([row1,row2],:)=A_([row2,row1],:);
        P_([row1,row2],:)=P_([row2,row1],:);
        L_([row1,row2],1:row1-1)=L_([row2,row1],1:row1-1);
    else
        % operations to process row2 with row1
        times=A_(row2,row1)/A_(row1,row1);
        A_(row2,:)=A_(row2,:)-times*A_(row1,:);
        L_(row2,row1)=times;
    end
end

% To check whether PA equals to LDU
function [status]=check_LDU(P,A,L,D,U)
% Inputs:
%   P: The returned permuatation matrix after decomposition
%   A: The original matrix for decomposition
%   L: The returned lower diagonal matrix after decomposition
%   D: The returned diagonal matrix after decomposition
%   U: The returned upper matrix after decomposition
% Outputs:
%   status: True when PA=LDU, else False
    status=false;
    Q1=P*A;
    Q2=L*D*U;
    if max(abs(Q1(:) - Q2(:))) < 0.00001
        status=true;
    end
end

function [U,S,V]=SVD_Decompose(A)
% Inputs:
%   A: The input m*n matrix that need SVD decomposition.
% Outputs:
%   U: The m*m matrix with columns are eigenvectors of AA'
%   S: The m*n matrix with the first r diagonal entries are square roots of
%   eigenvalues of AA'
%   V: The n*n matrix with columns are eigenvalues of A'A
%      The returned values make A=USV'
    [m,n]=size(A);
    S=zeros(m,n);
    [U,D]=eig(A*A');
    % sort the eigenvalues and eigenvectors to make the diagonal entries in
    % descending order
    [~,index]=sort(diag(D),'descend');
    U=U(:,index);
    [V,D]=eig(A'*A);
    [D_sort,index]=sort(diag(D),'descend');
    D=diag(D_sort);
    V=V(:,index);
    r=rank(D);
    S(1:r,1:r)=sqrt(D(1:r,1:r));
    V(:,1:r)=A'*U(:,1:r)/S(1:r,1:r);
end

% To check whether A equals to USV
function [status]=check_SVD(A,U,S,V)
% Inputs:
%   A: The original matrix for decomposition
%   U: The returned m*m matrix after SVD decompostion
%   S: The returned m*n diagonal matrix after SVD decomposition
%   V: The returned n*n matrix after SVD decomposition
% Outputs:
%   status: True when A=USV', else False
    status=false;
    Q=U*S*V';
    if max(abs(A(:) - Q(:))) < 0.00001
        status=true;
    end
end