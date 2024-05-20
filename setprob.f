      subroutine setprob
      implicit double precision (a-h,o-z)
      character*25 fname
	 
      common /cparam/ Re
c
c     # Set the material parameters for the INSE equations
c
      iunit = 7
      fname = 'setprob.data'

      call opendatafile(iunit, fname)
                
c     # Reynolds number
c      read(7,*) Re
	  
	  
      return
      end
