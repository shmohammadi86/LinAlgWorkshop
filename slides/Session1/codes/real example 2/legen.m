%%%% Compute the value of L_N^{(m)}(x)

function r=legen(N,m,x)


  for j=1:m+1

   if j==1

     s(1,j)=1;  s(2,j)=x;

      for k=1:N-1
         s(k+2,j)=((2*k+1)*x*s(k+1,j)-k*s(k,j))/(k+1);
      end

   else

     s(1,j)=0; 
     if j==2
        s(2,j)=1;
     else
        s(2,j)=0;
     end
    
     for k=1:N-1
        s(k+2,j)=(2*k+1)*s(k+1,j-1)+s(k,j);
     end
   end

  end


      r=s(N+1,m+1);

