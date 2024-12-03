
function [hatrho_J,hateta_J,var_source_Gamma1,var_source_Gamma2,jn,var_source_Gamma1_th,bound1,bound2,g31iJ,g31iJ_2] = eta_estimation2(e,omega,alpha,r,seed,eta_num)
    rng(seed);
    n=size(e,1);
    jn=floor(n^(alpha));
    Jn=zeros(jn,3);
    hatU1star_J=0;
    hatrho_J=0;
    hatU3star_J=0;
    doubles =[1,2;1,3;2,3];%nchoosek(1:r, 2);
    a3iJstar=zeros(n, 1);
    a1iJstar=zeros(n, 1);
    size_i_in_J=zeros(n,1);

    X3X1rho_ijk=zeros(jn,3);
    for l=1:jn
        
        ijk = sort(randperm(n, r),2); % Sample r unique integers from 1 to n
        Jn(l,:)=ijk;
        i=ijk(1);
        j=ijk(2);  
        k=ijk(3);
        size_i_in_J(i)=size_i_in_J(i)+1;
        size_i_in_J(j)=size_i_in_J(j)+1;
        size_i_in_J(k)=size_i_in_J(k)+1;
        if eta_num==3
            %hatU3star_J
            h3_star_ijk=(e(i,j).*e(i,k)+e(j,i).*e(j,k)+e(k,i).*e(k,j))/3;
        elseif eta_num==4
            h3_star_ijk=(e(i,j).*e(k,j)+e(j,i).*e(k,i)+e(i,k).*e(j,k))/3;
        elseif eta_num==5
            %hatU3star_J
            h3_star_ijk=(e(i,j).*e(j,k)+e(i,k).*e(k,j)+e(j,i).*e(i,k)+e(k,j).*e(j,i)+e(j,k).*e(k,i)+e(k,i).*e(i,j))./6;
        end
       
        %hatU1star_J and hatrho_J
        for t = 1:size(doubles, 1)
            choose_1=doubles(t,1);
            choose_2=doubles(t,2);
            i_1=ijk(choose_1); 
            i_2=ijk(choose_2);
            h1_star_i1i2=(e(i_1,i_2)+e(i_2,i_1))/2;
            if eta_num==2
                h2_star_i112=(e(i_1,i_2)*e(i_2,i_1));
                hatU3star_J=hatU3star_J+h2_star_i112/3/jn;
            end
            rho_i1i2=(omega(i_1,i_2)+omega(i_2,i_1))/2;
            hatU1star_J=hatU1star_J+h1_star_i1i2/nchoosek(r,2)/jn;
            hatrho_J=hatrho_J+rho_i1i2/nchoosek(r,2)/jn;
            %for reduced sample variance
            X3X1rho_ijk(l,2)=X3X1rho_ijk(l,2)+h1_star_i1i2/3;
            X3X1rho_ijk(l,3)=X3X1rho_ijk(l,3)+rho_i1i2/3;
            if eta_num==2
                X3X1rho_ijk(l,1)=X3X1rho_ijk(l,1)+h2_star_i112/3;
            end
            %g11iJstar*|I_3^{(l)}\in J:i\in I_3^{(l)}|
            if i_1==i
                a1iJstar(i)=a1iJstar(i)+h1_star_i1i2/2;
                if eta_num==2
                    a3iJstar(i)=a3iJstar(i)+h2_star_i112/2;
                end
            end
            if (i_1==j) | (i_2==j)
                a1iJstar(j)=a1iJstar(j)+h1_star_i1i2/2;
                if eta_num==2
                    a3iJstar(j)=a3iJstar(j)+h2_star_i112/2;
                end
            end
            if i_2==k
                a1iJstar(k)=a1iJstar(k)+h1_star_i1i2/2;
                if eta_num==2
                    a3iJstar(k)=a3iJstar(k)+h2_star_i112/2;
                end
            end   
        end
        
        %a31iJstar*|I_3^{(l)}\in J:i\in I_3^{(l)}|
        if (eta_num==3) | (eta_num==5)  | (eta_num==4)
            a3iJstar(i)=a3iJstar(i)+h3_star_ijk;
            a3iJstar(j)=a3iJstar(j)+h3_star_ijk;
            a3iJstar(k)=a3iJstar(k)+h3_star_ijk;
            %for reduced sample variance
            X3X1rho_ijk(l,1)=X3X1rho_ijk(l,1)+h3_star_ijk;
            %hatU3star_J
            hatU3star_J=hatU3star_J+h3_star_ijk/jn;
        end

    end
    safe_size_i_in_J = size_i_in_J;
    safe_size_i_in_J(safe_size_i_in_J == 0) = 1;
    hateta_J=max(0,hatrho_J^(-2)).*(hatU3star_J-hatU1star_J^2);
    g31iJstar=a3iJstar./safe_size_i_in_J - hatU3star_J;
    g11iJstar=a1iJstar./safe_size_i_in_J - hatU1star_J;
    g31iJstar(size_i_in_J == 0) = -hatU3star_J;
    g11iJstar(size_i_in_J == 0) = -hatU1star_J;


    %var_source_Gamma1_OLD =  max(0,hatrho_J^(-4)).*mean((3.*(g31-hatU3star_J)-4.*hatU1star_J.*(g11-hatU1star_J)).^2)/n;
    rval=3;
    if eta_num==3
        rval=3;
    elseif eta_num==5
        rval=3;
    elseif eta_num==2
        rval=2;
    end
    var_source_Gamma1 =  max(0,hatrho_J^(-4)).*mean((rval.*g31iJstar-4.*hatU1star_J.*g11iJstar).^2)/n;
    hatX3_ijk =  max(0,hatrho_J^(-2)).*(X3X1rho_ijk(:,1)-hatU3star_J-hatU1star_J.*X3X1rho_ijk(:,2)+hatU1star_J.^2);
    hatXrho3_ijk=2*max(0,hatrho_J^(-2)).*hateta_J.*(X3X1rho_ijk(:,3)-hatrho_J);
    var_source_Gamma2 = mean((hatX3_ijk-hatXrho3_ijk).^2)/jn;

    if var_source_Gamma1>n^(-3/2)*sqrt(log(n))+max(0,hatrho_J^(-2))*n^(-1/2-alpha)*(log(n))^(1/2)
        var_source_Gamma1_th=var_source_Gamma1;
    else
        var_source_Gamma1_th=0;
    end
    g31iJ= max(0,hatrho_J^(-2)).*(rval.*g31iJstar);
    g31iJ_2= max(0,hatrho_J^(-2)).*(rval.*g31iJstar-4.*hatU1star_J.*g11iJstar);

    bound1=n^(-3/2)*sqrt(log(n));
    bound2=max(0,hatrho_J^(-2))*n^(-1/2-alpha)*(log(n))^(1/2);
    %rate_Gamma2=rho^(-2)*n^(-alpha);
    %rate_gamma1_deg=(n^(-3/2)*sqrt(log(n))+rho^(-2)*n^(-1/2-alpha)*(log(n))^(1/2)),
    %rate_gamma1_nondeg=xi*n^(-1),
    %c=round(c,2),
    %rate_deg=rho^2*n^(alpha-3/2))

end







