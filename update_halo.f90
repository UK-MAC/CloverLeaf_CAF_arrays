!Crown Copyright 2012 AWE.
!
! This file is part of CloverLeaf.
!
! CloverLeaf is free software: you can redistribute it and/or modify it under 
! the terms of the GNU General Public License as published by the 
! Free Software Foundation, either version 3 of the License, or (at your option) 
! any later version.
!
! CloverLeaf is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details.
!
! You should have received a copy of the GNU General Public License along with 
! CloverLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief Driver for the halo updates
!>  @author Wayne Gaudin
!>  @details Invokes the kernels for the internal and external halo cells for
!>  the fields specified.

MODULE update_halo_module

CONTAINS

SUBROUTINE update_halo(fields,depth)

  USE clover_module
  USE update_halo_kernel_module

  IMPLICIT NONE

  INTEGER :: c,fields(NUM_FIELDS),depth

  CALL clover_exchange(fields,depth)

  DO c=1,chunks_per_task

    IF(chunks(c)%task.EQ.parallel%task) THEN

      IF(use_fortran_kernels)THEN
        CALL update_halo_kernel(chunks(c)%field%x_min,          &
                                chunks(c)%field%x_max,          &
                                chunks(c)%field%y_min,          &
                                chunks(c)%field%y_max,          &
                                chunks(c)%chunk_neighbours,     &
                                density0,       &
                                energy0,        &
                                pressure,       &
                                viscosity,      &
                                soundspeed,     &
                                density1,       &
                                energy1,        &
                                xvel0,          &
                                yvel0,          &
                                xvel1,          &
                                yvel1,          &
                                vol_flux_x,     &
                                vol_flux_y,     &
                                mass_flux_x,    &
                                mass_flux_y,    &
                                fields,                         &
                                depth                           )
      ELSEIF(use_C_kernels)THEN
        CALL update_halo_kernel_c(chunks(c)%field%x_min,        &
                                chunks(c)%field%x_max,          &
                                chunks(c)%field%y_min,          &
                                chunks(c)%field%y_max,          &
                                chunks(c)%chunk_neighbours,     &
                                density0,       &
                                energy0,        &
                                pressure,       &
                                viscosity,      &
                                soundspeed,     &
                                density1,       &
                                energy1,        &
                                xvel0,          &
                                yvel0,          &
                                xvel1,          &
                                yvel1,          &
                                vol_flux_x,     &
                                vol_flux_y,     &
                                mass_flux_x,    &
                                mass_flux_y,    &
                                fields,                         &
                                depth                           )
      ENDIF
    ENDIF

  ENDDO

END SUBROUTINE update_halo

END MODULE update_halo_module
