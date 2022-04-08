
   ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   ! USER INPUT : HOW MANY SOURCES DO WE NEED?
   ! -----------------------------------------
   self%numberOf_sources = 5
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
   ! shear localised at the top U(z)
   ! -----------------------------------------
   sourceDef = exp(-(1._dp-z) / self%sourceParams('shearL'))
   !sourceDef = z
   call self%add_source( definition = sourceDef,&
                         sourceIndex = 3)   ! this corresponds to source03
   ! -----------------------------------------
   ! shear localised at the top U'(z) (derivative)
   ! -----------------------------------------
   sourceDef = exp(-(1._dp-z) / self%sourceParams('shearL')) / &
                                  self%sourceParams('shearL')
   !sourceDef = 1._dp + 0._dp*z
   call self%add_source( definition = sourceDef,&
                         sourceIndex = 4)   ! this corresponds to source03
   ! -----------------------------------------
   ! shear localised at the top U'(z) (derivative)
   ! -----------------------------------------
   sourceDef = exp(-(1._dp-z) / self%sourceParams('shearL')) / &
                                  self%sourceParams('shearL')**2
   !sourceDef = 0._dp * z
   call self%add_source( definition = sourceDef,&
                         sourceIndex = 5)   ! this corresponds to source03
   ! -----------------------------------------
   ! /////////////////////////////////////////

   ! N.B.: both sources could have been defined as a unique source.
   ! To illustrate how to define multiple sources, I defined them separately.
