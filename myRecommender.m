function [ U, V ] = myRecommender( rateMatrix, lowRank )
    % Please type your name here:
    name = 'Gururani, Siddharth Kumar';
    disp(name); % Do not delete this line.

    % Parameters
    maxIter = 5000; % Choose your own.
    learningRate = 0.00058; % Choose your own.
    regularizer = 0.001; % Choose your own.
    
    % Random initialization:
    [n1, n2] = size(rateMatrix);
    U = rand(n1, lowRank) / lowRank;
    V = rand(n2, lowRank) / lowRank;
    
%     U_new = U;
%     V_new = V;
    % Gradient Descent:
    
    % IMPLEMENT YOUR CODE HERE.
    %new_error = sum(sum((rateMatrix - U*V').^2)) + regularizer*sum(sum(U.^2)) + regularizer*sum(sum(V.^2));
    %new_error = norm((U*V' - rateMatrix) .* (rateMatrix > 0), 'fro') / sqrt(nnz(rateMatrix > 0));
    new_error = (norm((U*V' - rateMatrix) .* (rateMatrix > 0), 'fro'))^2 + regularizer*sum(sum(U.^2)) + regularizer*sum(sum(V.^2));
    for iter = 1:maxIter
        % recompute U matrix
        U_old = U;
        V_old = V;
        old_error = new_error;
        for i = 1:n1
            populated = rateMatrix(i,:) > 0;
            err = U(i,:)*V';
            sum_term = populated .* (rateMatrix(i,:) - err);
            unregularized_derivative = -2 * sum_term * V;
            derivative = unregularized_derivative + 2*regularizer*U(i,:);
%             U_new(i,:) = U(i,:) + learningRate*derivative;
            U(i,:) = U(i,:) - learningRate*derivative;
        end
%             populated = rateMatrix(i,:) > 0;
%             err = U(i,:)*V';
%             sum_term = populated .* (rateMatrix(i,:) - err);
%             unregularized_derivative = -2 * sum_term * V;
%             derivative = unregularized_derivative + 2*regularizer*U(i,:);
% %             U_new(i,:) = U(i,:) + learningRate*derivative;
%             U(i,:) = U(i,:) - learningRate*derivative;
        % recompute V matrix
        for i = 1:n2
            populated = rateMatrix(:,i) > 0;
            populated = populated';
            err = V(i,:)*U';
            sum_term = populated .* (rateMatrix(:,i)' - err);
            unregularized_derivative = -2 * sum_term * U;
            derivative = unregularized_derivative + 2*regularizer*V(i,:);
%             V_new(i,:) = V(i,:) + learningRate*derivative;
            V(i,:) = V(i,:) - learningRate*derivative;
        end
%         new_error = sum(sum((rateMatrix - U*V').^2)) + regularizer*sum(sum(U.^2)) + regularizer*sum(sum(V.^2));
%         new_error = norm((U*V' - rateMatrix) .* (rateMatrix > 0), 'fro') / sqrt(nnz(rateMatrix > 0));
        new_error = (norm((U*V' - rateMatrix) .* (rateMatrix > 0), 'fro'))^2 + regularizer*sum(sum(U.^2)) + regularizer*sum(sum(V.^2));
        if(old_error - new_error < 0.05)
            U = U_old;
            V = V_old;
            break;
        end
%         U = U_new;
%         V = V_new;


%% Unoptimized
%     for i = 1:n1
%         for j = 1:n2
%             if(rateMatrix(i,j) ~= 0)
%                 err = rateMatrix(i,j) - dot(U(i,:),V(j,:));
%                 U(i,:) = U(i,:) + learningRate*(2*err*V(j,:) - 2*regularizer*U(i,:));
%                 V(j,:) = V(j,:) + learningRate*(2*err*U(i,:) - 2*regularizer*V(j,:));
%             end
%         end
%     end
% %     new_error = sum(sum((rateMatrix - U*V').^2)) + regularizer*sum(sum(U.^2)) + regularizer*sum(sum(V.^2));
%     new_error = norm((U*V' - rateMatrix) .* (rateMatrix > 0), 'fro') / sqrt(nnz(rateMatrix > 0));
%     if(old_error - new_error < 0.001)
%         break;
%     end

    end
    iter
end