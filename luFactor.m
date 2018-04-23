function [L, U, P] = luFactor(A)
% This function uses LU factorization to solve a matrix.  The luFactor
% function does not solve the equation, more it returns values for the L
% and U matrix.  Note that this is for square matrices. 
%   Inputs:
%       -A: the coefficient matrix
%   Outputs: 
%       -L: lower triangle matrix
%       -U: upper triangle matrix 
%       -P: the pivot matrix

if nargin < 1 %verify that a matrix is input. 
    error ('Input of the matrix is required');
end 

if ~isnumeric(A) %Verify that the matrix only contains numbers.
        error ('Matrix must contain only numbers')
end 
    
    [nRow,nCol] = size(A); %Assign the size of matrix A as values in the
    %nRow and nCol variables.
    if nRow ~= nCol %Verify the matrix is a square.
        error ('Matrix must be square')
    end 

    function pivot()    
        %This is a seperate function inside the main luFactor function.  It
        %performs the pivoting of the rows in order to perform LU
        %factorization.  It first identifies the row with the largest
        %coefficient for each column and then pivots the other elements to
        %get the coefficients in the right order. This does not perform any
        %of the solving that is involved with LU factorization but by
        %having this function in here, the main luFactor function can call
        %on the pivot() function when it needs to perform pivoting and
        %then solve as it needs to, then call upon pivot() to pivot the 
        %next row.
        [c,I] = max(abs(A(n:end,n)));
        I = I + (n-1);
        temp = A(n,:);
        A(n,:) = A(I,:);
        A(I,:) = temp;
        
        %Flips the appropriate L rows in order to perform correct pivoting 
        %on the matrices. 
        temp = L(n,:);
        L(n,:) = L(I,:);
        L(I,:) = temp;
        
        %Creates the matrix E to store values of P. Once all the pivoting
        %is performed, E*P = P. E is the permuation matrix.
        E(n,:) = 0;
        E(n,I) = 1;
        E(I,:) = 0;
        E(I,n) = 1;        
    end

    P = diag(ones(nRow,1)); %Fills the P matrix with 1's on the diagonal 
    %and 0's elsewhere to start LU factorization.
    U = zeros(nRow); %Creates a U matrix with 0's.
    L = zeros(nRow); %Creates an L matrix with 0's.

    for n = 1:nRow-1
        currentPivot = A(n,n); %After running the A matrix throught pivot()
        %function, it saves the iteration as a new variable in order to 
        %store all of the values used to perform factorization.
        
        E = diag(ones(nRow,1)); %Establishes the E matrix in order to 
        %store pivot values. The E matrix is filled with 1's on the
        %diagonal and 0's elsewhere.  Its size is based on however large
        %the input A matrix is. 
        
        maxPivot = max(A(n+1:end,n)); 
        if abs(currentPivot) < eps %Using eps to maintain accuracy and 
                                    %keep floating-point accuracy.
           if abs(maxPivot) < eps %Using eps to maintain accuracy and 
                                    %keep floating-point accuracy.
              error ('Cannot perform LU factorization with this matrix')
                        %If the errors calculated are too great or the
                        %matrix cannot be pivoted to be solved, then the
                        %luFactor function ends. 
           else
              pivot();  %Runs the pivot() function from above to perform 
                %the necessary pivoting on the matrices. 
           end
        else 
          if abs(currentPivot) < abs(maxPivot)
             if abs(currentPivot-maxPivot) >= 0 %Forces the accuracy of
                 % the factorization to be 100%.
                pivot();  %Runs pivot() function again when the accuracy 
                %is too far from being perfect. 
             end
          end
        end
        
        P = P*E;  %Calculate the new value for P as P*E.  E is filled in
            %during the pivot() function and then multiplied by P. 
        
        for i = n + 1:nRow
            L(i,n) = A(i,n)/A(n,n); %Takes the leading coefficients that 
            %need to be solved by elimination and divides them. Ex.
            %A21/A11. The value of this fraction then is placed in the
            %appropriate place on the L matrix.
            A(i,n) = 0;
            for j = n+1:nRow
                A(i,j) = A(i,j)-L(i,n)*A(n,j); %Using elimination to solve 
                %for the values of U.  The value calculated by dividing
                %A21/A11 is then multiplied by the row A11, and the row A21
                %is subtracted from it. The resulting values then are
                %stored as new values in the same A matrix until the matrix
                %has been completely factored out, then they will be
                %assigned to the U matrix.
            end
        end       
    end
   
    L = L + diag(ones(nRow,1)); %Creates a matrix the same size as input 
        %matrix A, filled with 1's on the diagonal and 0's elsewhere.
    P = P'; %Transposes P after it was multiplied by E.
    U = A; %Fills the U matrix with the calculated values from A, since
        %U*L = P*A.

L %Don't suppress the output of the variables in order to get them to be 
%displayed as matrices for the end user. 
U %Displays the U matrix.
P %Displays the P matrix.
    
end
