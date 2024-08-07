function [hatrho_J,hateta5_J,var_source_Gamma1,var_source_Gamma2,var_source_Gamma3,jn,var_source_Gamma1_2,var_source_Gamma1_3] = eta5_estimator_full(e,omega,alpha,r)
    % initial
    n=size(e,1);
    Jn1=nchoosek(1:n,2);
    Jn=nchoosek(1:n,3);
    jn=size(Jn,1);
    
    %get theta
    rowSums = sum(e, 2); % Column vector of row sums
    replicatedRowSums = repmat(rowSums, 1, n);
    thetastar_3 = (replicatedRowSums - e) / (n - 2); % is a n times n matrix;
    thetastar_3(logical(eye(size(thetastar_3)))) = 0; % Setting diagonal elements to 0 or other values
    
    %get some estimators
    doubles = nchoosek(1:r, 2);
    hatU1star_J=0;
    hatrho_J=0;
  
    %for t = 1:size(doubles, 1)
    i=1;
    j=2;
    t=1;
    k=3;
        %i=doubles(t,1);
        %j=doubles(t,2);
        idx_ij = Jn1(:,i)+(Jn1(:,j)-1)*n;
		idx_ji = Jn1(:,j)+(Jn1(:,i)-1)*n;
        e_ij=e(idx_ij);
        e_ji=e(idx_ji);
        omega_ij=omega(idx_ij);
        omega_ji=omega(idx_ji);
        hatU1star_J=hatU1star_J+mean((e_ij+e_ji)./2);  %/nchoosek(r,2)
        hatrho_J=hatrho_J+mean((omega_ij+omega_ji)./2);%/nchoosek(r,2)
    %end
    triplets = nchoosek(1:r, 3);
    hatU5star_J=0;
    
    %for t = 1:size(triplets, 1)
        %i=triplets(t,1);
        %j=triplets(t,2);
        %k=triplets(t,3);
        idx_ij = Jn(:,i)+(Jn(:,j)-1)*n;
		idx_jk = Jn(:,j)+(Jn(:,k)-1)*n;
		idx_ki = Jn(:,k)+(Jn(:,i)-1)*n;
        idx_kj = Jn(:,k)+(Jn(:,j)-1)*n;
		idx_ji = Jn(:,j)+(Jn(:,i)-1)*n;
		idx_ik = Jn(:,i)+(Jn(:,k)-1)*n;
        e_ij=e(idx_ij);
        e_ik=e(idx_ik);
        e_ji=e(idx_ji);
        e_jk=e(idx_jk);
        e_ki=e(idx_ki);
        e_kj=e(idx_kj);
        hatU5star_J=hatU5star_J+mean((e_ij.*e_jk+e_ik.*e_kj+e_ji.*e_ik+e_kj.*e_ji+e_jk.*e_ki+e_ki.*e_ij)./6);%/nchoosek(r,3);  
    %end
    
    hateta5star_J=hatU5star_J-hatU1star_J^2;
    hateta5_J=max(0,hatrho_J^(-2))*hateta5star_J;
    
    
    %get variance
    hatg_eta51_J=[];
    hatg_eta51_J_1=[];
    hatg_eta51_J_2=[];
    g51=[];
    g11=[];
    for idx_i = 1:n
        %g51(idx_i)=calculate_g51_J(e,idx_i,n,Jn);
        %g11(idx_i)=calculate_g11_J(e,idx_i,n,Jn);
        %hatg_eta51_J(idx_i)=max(0,hatrho_J^(-2)).*(3.*(g51(idx_i)-hatU5star_J)-4.*hatU1star_J.*(g11(idx_i)-hatU1star_J));
        hatg_eta51_J_2(idx_i)=max(0,hatrho_J^(-2)).*(3.*(calculate_g51_J(e,idx_i,n,nchoosek(1:n,3))-hatU5star_J)-4.*hatU1star_J.*(calculate_g11_J(e,idx_i,n,nchoosek(1:n,2))-hatU1star_J));
    end
   % var_source_Gamma1 = mean(hatg_eta51_J.^2)/n;
    var_source_Gamma1_2 = mean(hatg_eta51_J_2.^2)/n;
    %var_source_Gamma1_3 = max(0,hatrho_J^(-4)).*mean((3.*(g51-mean(g51))-4.*hatU1star_J.*(g11-mean(g11))).^2)/n;
    var_source_Gamma1=0;
    var_source_Gamma1_3=0;
    
     %get thetastar'
     
     Gamma2=0;
     var_source_Gamma2=0;
     if alpha>2
        for t = 1:size(doubles, 1)
            i=doubles(t,1);
            j=doubles(t,2);
            idx_ij = Jn(:,i)+(Jn(:,j)-1)*n;
            idx_ji = Jn(:,j)+(Jn(:,i)-1)*n;
            thetaprime_3_ij=thetastar_3(idx_ij).*max(0,hatrho_J ^(-2))-2.* max(0,hatrho_J ^(-1)).*hateta5_J/2;
            thetaprime_3_ji=thetastar_3(idx_ji).*max(0,hatrho_J ^(-2))-2.* max(0,hatrho_J ^(-1)).*hateta5_J/2;
            thetaprime_3_ij_mean=thetaprime_3_ij/jn/nchoosek(r,2);
            thetaprime_3_ji_mean=thetaprime_3_ji/jn/nchoosek(r,2);
        
            e_ij=e(idx_ij);
            e_ji=e(idx_ji);
            Gamma2=Gamma2+mean((thetaprime_3_ij.*e_ij+thetaprime_3_ji.*e_ji))/nchoosek(r,2);  
            var_source_Gamma2= var_source_Gamma2+sum(((thetaprime_3_ij_mean.^2).*e_ij+(thetaprime_3_ji_mean.^2).*e_ji));  
        end
     end

    
    var_source_Gamma3 = max(0,hatrho_J ^(-2*2)) * hatU5star_J / jn;
    

end
