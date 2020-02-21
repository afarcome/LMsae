      subroutine update_u4(Be,k,np1,MGa,np2,MDe,cN,cmax,U,M,X,mmax,nT,
     c                     Z,Y,n,TT,nnv,accU,tlrit,RR)

c definition of arguments
      integer k,np1,n,TT,np2,mmax,cmax,nT,cN(n,cmax),U(n,TT),M(n,TT)
      integer j,t,nnv(n)
      double precision Be(k,np1),accU,RR(k,n,TT)
      double precision MGa(k,np2),MDe(k,np2,k)
      double precision X(nT,np1-1),Y(nT)
      double precision Z(np2-k,n,TT),lrit(k),tlrit
c definition of interval variables
      integer j1,u1,m1,u0,u2,co,h1,rind(k)
      double precision lpv(mmax),llav(k),cova(np2),D(mmax,np1)
      double precision lav(k),fnv(k-1),lpiv(k),piv(k),rit(k),r,crit(k)
      double precision ml

c preliminaries 
      cova = 0; cova(1) = 1
      D(:,1) = 1
      tlrit = 0
      co = 0
      
c compute log-likelihoods
      do t = 1,TT
        do j = 1,n
          lrit = 0
c first time occasion
          if(t==1) then
            fnv = 0
            if(nnv(j)>0) then
              do h1=1,cmax
                j1 = cN(j,h1)
                if(j1>0) then
                  u2 = U(j1,1)-1
                  if(u2>0) fnv(u2) = fnv(u2)+1
                end if
              end do
              fnv = fnv/real(nnv(j))
            end if
            cova(2:(np2-k+1)) = Z(:,j,1)
            cova((np2-k+2):np2) = fnv
            llav = matmul(MGa,cova)
            lav = exp(llav)/sum(exp(llav))
            lrit = log(lav)
            fnv = 0
            if(nnv(j)>0) then
              do h1=1,cmax 
                j1 = cN(j,h1)
                if(j1>0) then
                  u2 = U(j1,2)-1
                  if(u2>0) fnv(u2) = fnv(u2)+1
                end if
              end do
              fnv = fnv/real(nnv(j))
            end if
            cova(2:(np2-k+1)) = Z(:,j,2)
            cova((np2-k+2):np2) = fnv
            do u1 = 1,k
              lpiv = matmul(MDe(:,:,u1),cova)
              piv = exp(lpiv)/sum(exp(lpiv))
              lrit(u1) = lrit(u1)+log(piv(U(j,2)))
            end do
c last time occasion
          else if(t==TT) then
            fnv = 0
            if(nnv(j)>0) then
              do h1=1,cmax
                j1 = cN(j,h1)
                if(j1>0) then
                  u2 = U(j1,TT)-1
                  if(u2>0) fnv(u2) = fnv(u2)+1
                end if
              end do
              fnv = fnv/real(nnv(j))
            end if
            cova(2:(np2-k+1)) = Z(:,j,TT)
            cova((np2-k+2):np2) = fnv
            lpiv = matmul(MDe(:,:,U(j,TT-1)),cova)
            piv = exp(lpiv)/sum(exp(lpiv))
            lrit = log(piv)
c other time occasions
          else
            fnv = 0
            if(nnv(j)>0) then
              do h1=1,cmax
                j1 = cN(j,h1)
                if(j1>0) then
                  u2 = U(j1,t)-1
                  if(u2>0) fnv(u2) = fnv(u2)+1
                end if
              end do
              fnv = fnv/real(nnv(j))
            end if
            cova(2:(np2-k+1)) = Z(:,j,t)
            cova((np2-k+2):np2) = fnv
            lpiv = matmul(MDe(:,:,U(j,t-1)),cova)
            piv = exp(lpiv)/sum(exp(lpiv))
            lrit = log(piv)
            fnv = 0
            if(nnv(j)>0) then
              do h1=1,cmax
                j1 = cN(j,h1)
                if(j1>0) then
                  u2 = U(j1,t+1)-1
                  if(u2>0) fnv(u2) = fnv(u2)+1
                end if
              end do
              fnv = fnv/real(nnv(j))
            end if
            cova(2:(np2-k+1)) = Z(:,j,t+1)
            cova((np2-k+2):np2) = fnv
            do u1 = 1,k
              lpiv = matmul(MDe(:,:,u1),cova)
              piv = exp(lpiv)/sum(exp(lpiv))
              lrit(u1) = lrit(u1)+log(piv(U(j,t+1)))
            end do
          end if
      
c contribution of response variables
          m1 = M(j,t)
          if(m1>0) then
            do u1 = 1,k
              D(1:m1,2:np1) = X((co+1):(co+m1),:)
              lpv(1:m1) = matmul(D(1:m1,:),Be(u1,:))
              lrit(u1) = lrit(u1)+sum(Y((co+1):(co+m1))*lpv(1:m1))
     c                   -sum(log(1+exp(lpv(1:m1))))
            end do
            co = co+M(j,t)
          end if

c update U
          ml = maxval(lrit)
          rit = exp(lrit-ml)/sum(exp(lrit-ml))
          RR(:,j,t) = lrit
          r = rand(0)
          crit = rit
          rind = 0
          if(r>crit(1)) rind(1) = 1
          do u1 = 2,k
            crit(u1) = crit(u1-1)+rit(u1)
            if(r>crit(u1)) rind(u1) = 1
          end do
          u0 = U(j,t)
          U(j,t) = 1+sum(rind)
          if(u0.ne.U(j,t)) accU = accU+1/real(n*TT)
          tlrit = tlrit + log(rit(U(j,t)))
          
        end do
      end do
      
      	
      end