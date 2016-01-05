function [ min_L, min_R ] = myRecommender_grid( rateMatrix, lowRank )
    % Please type your name here:
    name = 'Gururani, Siddharth Kumar';
    disp(name); % Do not delete this line.

    % Parameters
    maxIter = 1000; % Choose your own.
%     learningRate = 0.0002; % Choose your own.
%     regularizer = 0.002; % Choose your own.
    learningRate_set = logspace(-3,-5);
    regularizer_set = logspace(-2,2);
    min_er = Inf;
    for a = 3:50
        for b = 1:50
            learningRate = learningRate_set(a);
            regularizer = regularizer_set(b);
    % Random initialization:
    [n1, n2] = size(rateMatrix);
    U = rand(n1, lowRank) / lowRank;
    V = rand(n2, lowRank) / lowRank;
    
    U_old = U;
    V_old = V;
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
        if(old_error - new_error < 0.1)
            U = U_old;
            V = V_old;
            break;
        end
%         U = U_new;
%         V = V_new;



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
    logtime = toc;
    if(min_er > new_error)
        min_L = learningRate;
        min_R = regularizer;
    end
    trainRMSE = norm((U*V' - rateMatrix) .* (rateMatrix > 0), 'fro') / sqrt(nnz(rateMatrix > 0));
    fprintf('SVD-%d\t%.4f\t%f\t%f min_L %f\t min_R %f\n', 1, trainRMSE, iter, logtime, min_L, min_R);
        end
    end
end