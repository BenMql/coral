! USER: specify how many masks should be created:
N_masks = 2


! do not modify: 
if (allocated(penalisation_mask)) deAllocate (penalisation_mask)
allocate ( penalisation_mask (N_masks))
pi = 3.14159265358979311599796346854 ! in case you need it...

! USER: you now have to define as many masks as N_masks defined above.
! For each mask, increment i_mask and define the penalisation function
!
! exemple 1: penalisation function for domain above a given surface
!             z >  - (0.9 + 0.1*cos( 5*(x+y)*2*pi / lx))
!
 i_mask = i_mask + 1 
 penalisation_mask(i_mask) = z - (0.9 + 0.1*cos( 5*(x+y)*2*pi / lx))

!
! exemple 2: penalisation function for domain below a given surface:
!             z >  - (0.9 + 0.1*cos( 5*(x+y)*2*pi / lx))
!
 i_mask = i_mask + 1 
 penalisation_mask(i_mask) =  (0.4_dp * exp(                  &
                          - ((6._dp + cos(2*pi*x/lx) )**2)    &
                                           )                  &
                                       / exp( -25._dp)) -   z 
