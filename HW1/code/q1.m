clear;

% Test cases
A=[[1,-2,1];[1,2,2];[2,3,4]];
%A=[[1,1,0];[1,1,2];[4,2,3]];
%A=[[2,2,5];[1,1,5];[3,2,5]];
%A=[[-1,1,0,0];[-1,0,1,0];[0,-4,1,0];[0,-1,0,1]];
%A=[[0,-4,1,0];[-1,0,1,0];[-1,1,0,0];[0,-1,0,1]];

% Decompose A with PA=LDU
[P,L,D,U]=LDU_Decompose(A);

% Check the correctness of the deposition
status=check(P,A,L,D,U);


function [P,L,D,U]=LDU_Decompose(A)
% Inputs: 
%   A: The input matrix that is square and invertible
%
% Outputs:
%   P: The permutation matrix that is used for shifting the rows of A
%   L: The decomposition result L, which is a lower triangular matrix
%   D: The decomposition result D, which is a diagonal matrix
%   U: The decomposition result U, which is a upper triangular matrix
%      The returned matrixes make PA = LDU
    [m,~]=size(A);
    A_=A;
    P=eye(m);
    L=eye(m);
    for k=1:m
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
    D=diag(diag(A_));
    U=D\A_;
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
function [status]=check(P,A,L,D,U)
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