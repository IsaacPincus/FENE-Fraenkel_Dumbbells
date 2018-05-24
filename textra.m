clear variables

filename = 'to_textra.dat';
data = tdfread(filename);
data = cell2mat(struct2cell(data));

limit = 0.25;

x = data(:,1);
y = data(:,2);
dy = data(:,3);

iflag = 1;
min_err = 100;
intercept = 100;

final_coeff = [];
final_err = [];
pfin = [];

for i = 1:length(x)-1
    xi = x(i:end);
    yi = y(i:end);
    dyi = dy(i:end);
    
    for j = 1:(length(xi)-1)
        if iflag == 1
            p = j:-1:0;
        elseif iflag == 0
            p = [j:-1:2,0];
        end
        [coeff, error, chi2val] = pfit2(p, xi, yi, dyi);
        storage(i, j, 1) = coeff(end);
        storage(i, j, 2) = error(end);
        storage(i, j, 3) = chi2val;
        storage(i, j, 4) = chi2inv(1-limit, length(xi) - length(p));
        
        if chi2val < storage(i,j,4)
            if error(end) < min_err
                min_err = error(end);
                
                final_coeff = coeff;
                final_err = error;
                pfin = p;
            end
        end
    end
end

figure();
hold on
plot([0;x], bsxfun(@power, [0;x], pfin)*final_coeff, 'b-');
plot(x, y, 'ro');
title('TEXTRA minimum error extrapolation')
xlabel('\Deltat, timestep width')
ylabel('Extrapolated quantity')
hold off

function [coeff, error, Q] = pfit2(p, x, y, dy)
    V = bsxfun(@power, x, p);

    A = V./dy;
    b = y./dy;
    alpha = A'*A;
    beta = A'*b;
    coeff = alpha^-1*beta;
    error = sqrt(diag(alpha^-1));
    
    pi = V*coeff;

    Q = sum((y-pi).^2./dy.^2);
end

% function [coeff, error, passes_test] = pfit(p, x, y, dy, limit)
%     V = bsxfun(@power, x, p);
% 
%     [coeff, error] = lscov(V./dy, y./dy);
% 
%     pi = V*coeff;
% 
%     Q = sum((y-pi).^2./dy.^2);
%     chi2 = chi2inv(1-limit, length(x) - length(p));
%     prob = 1 - chi2cdf(Q, length(x) - length(p));
%     passes_test = (Q < chi2);
%     
% end