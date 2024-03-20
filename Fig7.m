% Bypassing the quadrature exactness assumption of hyperinterpolation on the sphere
% by C. An and H.-N. Wu
% Journal of Complexity (80) 2024, 101789

% MATLAB Codes written by H.-N. Wu in 2022

% Please add the sphere_approx_toolbox_v3.0 onto path

clear
close all

point_idx = 1;
Ls = [2 4 6 8 10 12 14 16 18];
Ks = 2;
mtps = [1 2 3 4 5];
i_f = 2;

% Validation point set
Xt = get_Xt( );

for ell = 1:max(size(mtps));
    multi = mtps(ell);

    for i = 1:max(size(Ls))
        L = Ls(i);
        m =multi* ceil((L+1)^2 * L^(2+2/(i_f+3/2)));
        fprintf(['(beta,n) = (',num2str(multi),', ',num2str(L),')\n'])

        % function to be approximated
        funtxt_ce = {'nrWend_k0','nrWend_k1','nrWend_k2','nrWend_k3','nrWend_k4'};
        funtxt = funtxt_ce{i_f+1};
        switch funtxt
            case 'nrWend_k0'
                rbf_k = 0;
            case 'nrWend_k1'
                rbf_k = 1;
            case 'nrWend_k2'
                rbf_k = 2;
            case 'nrWend_k3'
                rbf_k = 3;
            case 'nrWend_k4'
                rbf_k = 4;
        end
        switch funtxt
            case {'nrWend_k0','nrWend_k1','nrWend_k2','nrWend_k3','nrWend_k4'}
                func = @rbf_nr;
        end
        k = rbf_k;
        delta = (3*k+3)*gamma(k+1/2)/(2*gamma(k+1));
        ft = func(Xt',rbf_k);


        switch point_idx
            case 1 % Equal area points
                X_k = sca_data_cap(2, 1 ,m, 0);
                X_k = X_k';
            case 2 % Minimal energy points
                X_k = loadME( t, (t+1)^2 );
            case 3
                X_k = loadMD( t, (t+1)^2 );
            case 4
                X_k = loadStd(t,(t+1)^2);
        end

        % function sampling
        f = func(X_k,rbf_k);
        Y_L = getQ( X_k, L )';
        alpha = 4*pi*Y_L*f/m;

        % approximation polynomials
        L = sqrt(length(alpha))-1;
        Yt_eqp = get_Yt( L, Xt );
        pt_eqphyper = Yt_eqp * alpha;
        error(i,ell) = sqrt(4*pi/50000*sum((ft-pt_eqphyper).^2));
    end
end

semilogy(Ls,error(:,1),'marker','square', 'LineWidth',2,...
    'Color',[0 0.4470 0.7410],...
    'MarkerSize',10,...
    'MarkerEdgeColor',[0 0.4470 0.7410],...
    'MarkerFaceColor',[0 0.4470 0.7410])
hold on
semilogy(Ls,error(:,2),'marker','pentagram', 'LineWidth',2,...
    'Color',[0.8500 0.3250 0.0980]	,...
    'MarkerSize',10,...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980]	,...
    'MarkerFaceColor',[0.8500 0.3250 0.0980]	)
semilogy(Ls,error(:,3),'marker','o', 'LineWidth',2,...
    'Color',[0.9290 0.6940 0.1250]		,...
    'MarkerSize',10,...
    'MarkerEdgeColor',[0.9290 0.6940 0.1250]	,...
    'MarkerFaceColor',[0.9290 0.6940 0.1250]	)
semilogy(Ls,error(:,4),'marker','d', 'LineWidth',2,...
    'Color',[0.4660 0.6740 0.1880]		,...
    'MarkerSize',10,...
    'MarkerEdgeColor',[0.4660 0.6740 0.1880]			,...
    'MarkerFaceColor',[0.4660 0.6740 0.1880]		)
semilogy(Ls,error(:,5),'marker','^', 'LineWidth',2,...
    'Color',[0.3010 0.7450 0.9330]			,...
    'MarkerSize',10,...
    'MarkerEdgeColor',[0.3010 0.7450 0.9330]				,...
    'MarkerFaceColor',[0.3010 0.7450 0.9330]			)
hold on
semilogy(Ls,Ls.^(-(i_f+3/2)),'--k', 'LineWidth',3)

set(gca, 'fontsize', 24)
xlabel('degree $n$','interpreter','latex','fontsize',33), ylabel('error','interpreter','latex','fontsize',36)

set(legend( '$\beta = 1 $','$\beta = 2$', '$\beta = 3 $', '$\beta = 4 $', '$\beta = 5 $','$n^{-(\sigma+3/2)}$'),'interpreter','latex','fontsize',24)
title('$m=\beta\lceil(n+1)^2n^{2+\frac{2}{\sigma+3/2}}\rceil$ with $\sigma=2$','interpreter','latex','fontsize',30)
box on, grid on,


