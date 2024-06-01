function Jn=getsubset(n,alpha,r)
    jn=floor(n^(alpha));
    Jn=[];
    count_valid=0;
	while(count_valid<jn)
		temp_mat = randi([1,n], floor(jn/10), r);
		unique_columns = logical(ones(floor(jn/10),1));
		for(rr1 = 1:(r-1))
			for(rr2 = (rr1+1):r)
				unique_columns = unique_columns & (temp_mat(:,rr1)~=temp_mat(:,rr2));
			end
		end
		count_valid = count_valid + sum(unique_columns);
		Jn = [Jn; temp_mat(unique_columns,:)];
	end
	Jn = Jn(1:jn,:);
    Jn=sort(Jn,2);
end