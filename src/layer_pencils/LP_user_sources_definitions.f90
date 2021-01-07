
   ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   ! USER INPUT : HOW MANY SOURCES DO WE NEED?
   ! -----------------------------------------
   self%numberOf_sources = 2
   ! -----------------------------------------
   ! /////////////////////////////////////////
   
   call self%allocate_sources() ! do not modify or displace
   
   ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   ! USER INPUT : SOURCES DEFINITIONS            
   ! -----------------------------------------
   ! Heating localised at the bottom 
   ! -----------------------------------------
   sourceDef = self%sourceParams('HovEll') * exp(-z * self%sourceParams('HovEll') ) / &
                 (1._dp - exp(-self%sourceParams('HovEll') ) ) - 1._dp
   call self%add_source( definition = sourceDef,&
                         sourceIndex = 1)   ! this corresponds to source01
   ! -----------------------------------------
   ! Homogeneous cooling (more concise...)
   ! -----------------------------------------
   sourceDef = -1._dp 
   call self%add_source( definition = sourceDef,&
                         sourceIndex = 2)   ! this corresponds to source02
   ! -----------------------------------------
   ! /////////////////////////////////////////

   ! N.B.: both sources could have been defined as a unique source.
   ! To illustrate how to define multiple sources, I defined them separately.
