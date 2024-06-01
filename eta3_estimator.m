function [hatrho_J,hateta3_J,hateta3_full,var_source_Gamma1,var_source_gamma1_2_2,var_source_Gamma3,jn,var_source_Gamma1_2,var_source_Gamma1_full] = eta3_estimator(e,omega,alpha,r,Jn)
    % initial
    n=size(e,1);
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
    for t = 1:size(doubles, 1)
        i=doubles(t,1);
        j=doubles(t,2);
        idx_ij = Jn(:,i)+(Jn(:,j)-1)*n;
		idx_ji = Jn(:,j)+(Jn(:,i)-1)*n;
        e_ij=e(idx_ij);
        e_ji=e(idx_ji);
        omega_ij=omega(idx_ij);
        omega_ji=omega(idx_ji);
        hatU1star_J=hatU1star_J+mean((e_ij+e_ji)./2/nchoosek(r,2));  %
        hatrho_J=hatrho_J+mean((omega_ij+omega_ji)./2/nchoosek(r,2));%
    end
    triplets = nchoosek(1:r, 3);
    hatU3star_J=0;
    var=0;
    for t = 1:size(triplets, 1)
        i=triplets(t,1);
        j=triplets(t,2);
        k=triplets(t,3);
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
        hatU3star_J=hatU3star_J+mean((e_ij.*e_ik+e_ji.*e_jk+e_ki.*e_kj)./3/nchoosek(r,3));%;  
        %var=var+mean( ((e_ij.*e_ik+e_ji.*e_jk+e_ki.*e_kj)./3-hatU3star_J).^2)/jn/(nchoosek(r,3)^2);

    end

    %for t = 1:size(triplets, 1)
    %    i=triplets(t,1);
    %    j=triplets(t,2);
    %    k=triplets(t,3);
    %    idx_ij = Jn(:,i)+(Jn(:,j)-1)*n;
	%	idx_jk = Jn(:,j)+(Jn(:,k)-1)*n;
	%	idx_ki = Jn(:,k)+(Jn(:,i)-1)*n;
    %   idx_kj = Jn(:,k)+(Jn(:,j)-1)*n;
	%	idx_ji = Jn(:,j)+(Jn(:,i)-1)*n;
	%	idx_ik = Jn(:,i)+(Jn(:,k)-1)*n;
    %    e_ij=e(idx_ij);
    %    e_ik=e(idx_ik);
    %   e_ji=e(idx_ji);
    %    e_jk=e(idx_jk);
    %    e_ki=e(idx_ki);
    %    e_kj=e(idx_kj);
    %    %hatU3star_J=hatU3star_J+mean((e_ij.*e_ik+e_ji.*e_jk+e_ki.*e_kj)./3/nchoosek(r,3));%;  
    %    var=var+mean( ((e_ij.*e_ik+e_ji.*e_jk+e_ki.*e_kj)./3-hatU3star_J).^2)/jn/(nchoosek(r,3)^2);

    %end
    
    hateta3star_J=hatU3star_J-hatU1star_J^2;
    hateta3_J=max(0,hatrho_J^(-2))*hateta3star_J;
    
    sizealpha=floor(n^(alpha-1));
    var3=[];
    var=[];
    var1=[];
    tildemeansquare3=[];
    tildemeansquare1=[];
    tildemeansquare=[];
    id=0;
    e(logical(eye(size(e)))) = 0;
    for i=0:(n-1)
        for d=1:sizealpha
            id=id+1;
            idx_1=i:d:(i+(r-1)*d);
            idx_1=mod(idx_1,n)+1;
            idx_2=i:(-d):(i-(r-1)*d);
            idx_2=mod(idx_2,n)+1;
            idx_3=(i+r*d):(d):(i+(2*r-1)*d);
            idx_3=mod(idx_3,n)+1;
            idx1_ij=idx_1(1)+(idx_1(2)-1)*n;
            idx1_ji= idx_1(2)+(idx_1(1)-1)*n;
            idx1_ik=idx_1(1)+(idx_1(3)-1)*n;
            idx1_ki= idx_1(3)+(idx_1(1)-1)*n;
            idx1_jk=idx_1(2)+(idx_1(3)-1)*n;
            idx1_kj=idx_1(3)+(idx_1(2)-1)*n;
            idx2_ij=idx_2(1)+(idx_2(2)-1)*n;
            idx2_ji= idx_2(2)+(idx_2(1)-1)*n;
            idx2_ik=idx_2(1)+(idx_2(3)-1)*n;
            idx2_ki= idx_2(3)+(idx_2(1)-1)*n;
            idx2_jk=idx_2(2)+(idx_2(3)-1)*n;
            idx2_kj= idx_2(3)+(idx_2(2)-1)*n;
            idx3_ij=idx_3(1)+(idx_3(2)-1)*n;
            idx3_ji= idx_3(2)+(idx_3(1)-1)*n;
            idx3_ik=idx_3(1)+(idx_3(3)-1)*n;
            idx3_ki= idx_3(3)+(idx_3(1)-1)*n;
            idx3_jk=idx_3(2)+(idx_3(3)-1)*n;
            idx3_kj=idx_3(3)+(idx_3(2)-1)*n;
            e1_ij=e(idx1_ij);
            e1_ik=e(idx1_ik);
            e1_ji=e(idx1_ji);
            e1_jk=e(idx1_jk);
            e1_ki=e(idx1_ki);
            e1_kj=e(idx1_kj);
            e2_ij=e(idx2_ij);
            e2_ik=e(idx2_ik);
            e2_ji=e(idx2_ji);
            e2_jk=e(idx2_jk);
            e2_ki=e(idx2_ki);
            e2_kj=e(idx2_kj);
            e3_ij=e(idx3_ij);
            e3_ik=e(idx3_ik);
            e3_ji=e(idx3_ji);
            e3_jk=e(idx3_jk);
            e3_ki=e(idx3_ki);
            e3_kj=e(idx3_kj);
            varfirst=3*((e1_ij.*e1_ik+e1_ji.*e1_jk+e1_ki.*e1_kj)./3)-4*hatU1star_J*((e1_ij+e1_ji)/6+(e1_ik+e1_ki)/6+(e1_kj+e1_jk)/6);
            varsecond=3*((e2_ij.*e2_ik+e2_ji.*e2_jk+e2_ki.*e2_kj)./3)-4*hatU1star_J*((e2_ij+e2_ji)/6+(e2_ik+e2_ki)/6+(e2_kj+e2_jk)/6);
            varthird=3*((e3_ij.*e3_ik+e3_ji.*e3_jk+e3_ki.*e3_kj)./3)-4*hatU1star_J*((e3_ij+e3_ji)/6+(e3_ik+e3_ki)/6+(e3_kj+e3_jk)/6);
            tildemeansquare(id)=varfirst*varthird;
            %tildemeansquare3(id)=((e1_ij.*e1_ik+e1_ji.*e1_jk+e1_ki.*e1_kj)./3)*((e3_ij.*e3_ik+e3_ji.*e3_jk+e3_ki.*e3_kj)./3);
            var(id)=varfirst*varsecond;
            %tildemeansquare1(id)=((e1_ij+e1_ji)/6+(e1_ik+e1_ki)/6+(e1_kj+e1_jk)/6)*((e3_ij+e3_ji)/6+(e3_ik+e3_ki)/6+(e3_kj+e3_jk)/6);
            %tildemeansquare(id)=(9*tildemeansquare3(id)-16*hatU1star_J^2*tildemeansquare1(id))/jn;
        end
    end

    var_V_J_star=mean( ((e_ij.*e_ik+e_ji.*e_jk+e_ki.*e_kj)./3-hatU3star_J-hatU1star_J.*((e_ij+e_ji)./2-hatU1star_J)).^2)/jn;
    musquare=(sum(tildemeansquare));
    var_V_J_star_2=abs(sum(var-sum(tildemeansquare)/jn)/(jn^2)*sizealpha);

    
    
    %get variance
    g31=[];
    g11=[];
    g31_full=[];
    g11_full=[];
    for idx_i = 1:n
        g31(idx_i)=calculate_g31_J(e,idx_i,n,Jn);
        g11(idx_i)=calculate_g11_J(e,idx_i,n,Jn);
        g31_full(idx_i)=calculate_g31_J(e,idx_i,n,nchoosek(1:n,3));
        g11_full(idx_i)=calculate_g11_J(e,idx_i,n,nchoosek(1:n,2));
    end
    var_source_Gamma1 =  max(0,hatrho_J^(-4)).*mean((3.*(g31-hatU3star_J)-4.*hatU1star_J.*(g11-hatU1star_J)).^2)/n;
    var_source_Gamma1_2 =  max(0,hatrho_J^(-4)).*(var_V_J_star);
    var_source_Gamma1_full = max(0,hatrho_J^(-4)).*mean((3.*(g31_full-mean(g31_full))-4.*mean(g11_full).*(g11_full-mean(g11_full))).^2)/n;
    
    hateta3_full=mean(g31_full)-mean(g11_full)^2;
    
    var_source_gamma1_2_2=max(0,hatrho_J^(-4)).*var_V_J_star_2;
     %get thetastar'
     
     Gamma2=0;
     var_source_Gamma2=0;
     if alpha>2
        for t = 1:size(doubles, 1)
            i=doubles(t,1);
            j=doubles(t,2);
            idx_ij = Jn(:,i)+(Jn(:,j)-1)*n;
            idx_ji = Jn(:,j)+(Jn(:,i)-1)*n;
            thetaprime_3_ij=thetastar_3(idx_ij).*max(0,hatrho_J ^(-2))-2.* max(0,hatrho_J ^(-1)).*hateta3_J/2;
            thetaprime_3_ji=thetastar_3(idx_ji).*max(0,hatrho_J ^(-2))-2.* max(0,hatrho_J ^(-1)).*hateta3_J/2;
            thetaprime_3_ij_mean=thetaprime_3_ij/jn/nchoosek(r,2);
            thetaprime_3_ji_mean=thetaprime_3_ji/jn/nchoosek(r,2);
        
            e_ij=e(idx_ij);
            e_ji=e(idx_ji);
            Gamma2=Gamma2+mean((thetaprime_3_ij.*e_ij+thetaprime_3_ji.*e_ji))/nchoosek(r,2);  
            var_source_Gamma2= var_source_Gamma2+sum(((thetaprime_3_ij_mean.^2).*e_ij+(thetaprime_3_ji_mean.^2).*e_ji));  
        end
     end

    
    var_source_Gamma3 = max(0,hatrho_J ^(-2*2)) * hatU3star_J / jn;
    

end
