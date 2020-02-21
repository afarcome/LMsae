      subroutine lk_ga2(MGa,k,np2,MGas,cN,cmax,n,U,TT,Z,nnv,lk,lks)

c definition of arguments
      integer k,np2,TT,n,cmax,cN(n,cmax),U(n,TT),nnv(n)
      double precision MGa(k,np2),MGas(k,np2),Z(np2-k,n,TT),lk,lks
c definition of interval variables
      integer j,j1,u2,h1
      double precision fnv(k-1),cova(np2),llav(k),lav(k)

c preliminaries
      cova = 0; cova(1) = 1

c compute log-likelihoods
      lk = 0; lks = 0
      do j=1,n
        fnv = 0
        if(nnv(j)>0) then
          do h1=1,cmax
            j1 = cN(j,h1)
            if(j1>0) then
              u2 = U(j1,1)-1
              if(u2>0) fnv(u2) = fnv(u2)+1
            end if
          end do
          fnv = fnv/nnv(j)
        end if
        cova(2:(np2-k+1)) = Z(:,j,1)
        cova((np2-k+2):np2) = fnv
        llav = matmul(MGa,cova)
        lav = exp(llav)/sum(exp(llav))
        lk = lk+log(lav(U(j,1)))
        llav = matmul(MGas,cova)
        lav = exp(llav)/sum(exp(llav))
        lks = lks+log(lav(U(j,1)))
      end do

      end