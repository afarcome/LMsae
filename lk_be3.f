      subroutine lk_be3(Be,Bes,U,M,X,Y,n,TT,k,np1,nT,mmax,lkv,lksv)

      integer np1,n,TT,k,j,t,nT,mmax,U(n,TT),M(n,TT),u1,m1,co
      double precision X(nT,np1-1),Y(nT)
      double precision Be(k,np1),Bes(k,np1),lkv(k),lksv(k)
      double precision lpv(mmax),lpsv(mmax)
      double precision D(mmax,np1)

c compute log-likelihoods
      lkv = 0; lksv = 0
      D = 0; D(:,1) = 1
      co = 0
      do t=1,TT
        do j=1,n
          m1 = M(j,t)
          if(m1>0) then
            u1 = U(j,t)
            D(1:m1,2:np1) = X((co+1):(co+m1),:)
            lpv(1:m1) = matmul(D(1:m1,:),Be(u1,:))
            lkv(u1) = lkv(u1)+sum(Y((co+1):(co+m1))*lpv(1:m1))
     c                -sum(log(1+exp(lpv(1:m1))))
            lpsv(1:m1) = matmul(D(1:m1,:),Bes(u1,:))
            lksv(u1) = lksv(u1)+sum(Y((co+1):(co+m1))*lpsv(1:m1))
     c                 -sum(log(1+exp(lpsv(1:m1))))
            co = co+M(j,t)
          end if
        end do
      end do

      end