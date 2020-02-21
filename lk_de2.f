      subroutine lk_de2(MDe,k,np2,MDes,C,n,U,TT,Z,nnv,up,lk,lks)

c definition of arguments
      integer k,np2,TT,n,C(n,n),U(n,TT),nnv(n),up
      double precision MDe(k,np2,k),MDes(k,np2,k),Z(np2-k,n,TT),lk,lks
c definition of interval variables
      integer j,j1,t,u2
      double precision fnv(k-1),cova(np2),lpiv(k),piv(k)

c preliminaries
      cova = 0; cova(1) = 1

c compute log-likelihoods
      lk = 0; lks = 0
      do j=1,n
        do t=2,TT
          if(U(j,t-1)==up) then
            fnv = 0
            if(nnv(j)>0) then
              do j1=1,n
                u2 = U(j1,t)-1
                if(C(j,j1)==1.and.u2>0) fnv(u2) = fnv(u2)+1
              end do
              fnv = fnv/nnv(j)
            end if
            cova(2:(np2-k+1)) = Z(:,j,t)
            cova((np2-k+2):np2) = fnv
            lpiv = matmul(MDe(:,:,up),cova)
            piv = exp(lpiv)/sum(exp(lpiv))
            lk = lk+log(piv(U(j,t)))
            lpiv = matmul(MDes(:,:,up),cova)
            piv = exp(lpiv)/sum(exp(lpiv))
            lks = lks+log(piv(U(j,t)))
          end if
        end do
      end do

      end