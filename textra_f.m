function [Q, dQ] = textra_f(filename, iflag, limit)

    data = tdfread(filename);
    data = cell2mat(struct2cell(data));

    x = data(:,1);
    y = data(:,2);
    dy = data(:,3);

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
                if chi2val/(length(xi) - length(p)) < min_err
                    min_err = chi2val/(length(xi) - length(p));

                    final_coeff = coeff;
                    final_err = error;
                    pfin = p;
                end
            end
        end
    end
    
    Q = final_coeff(end);
    dQ = final_err(end);
end

function [coeff, error, Q] = pfit2(p, x, y, dy)
    %Implements polyfit scheme in 15.4 Numerical Methods 3rd Ed
        % p is horizontal array containing polynomial powers
        % x is timestep widths (or general x-data)
        % y is dependent data
        % dy is standard error in each y point
    
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