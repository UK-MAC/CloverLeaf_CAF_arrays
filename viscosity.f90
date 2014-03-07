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

!>  @brief Driver for the viscosity kernels
!>  @author Wayne Gaudin
!>  @details Selects the user specified kernel to caluclate the artificial 
!>  viscosity.

MODULE viscosity_module

CONTAINS

SUBROUTINE viscosity_driver()

  USE clover_module
  USE viscosity_kernel_module
  
  IMPLICIT NONE

  INTEGER :: c

  DO c=1,chunks_per_task

    IF(chunks(c)%task.EQ.parallel%task) THEN

      IF(use_fortran_kernels)THEN
        CALL viscosity_kernel(chunks(c)%field%x_min,                   &
                            chunks(c)%field%x_max,                     &
                            chunks(c)%field%y_min,                     &
                            chunks(c)%field%y_max,                     &
                            chunks(c)%field%celldx,                    &
                            chunks(c)%field%celldy,                    &
                            density0,                  &
                            pressure,                  &
                            viscosity,                 &
                            xvel0,                     &
                            yvel0                      )
      ELSEIF(use_C_kernels)THEN
        CALL viscosity_kernel_c(chunks(c)%field%x_min,                 &
                            chunks(c)%field%x_max,                     &
                            chunks(c)%field%y_min,                     &
                            chunks(c)%field%y_max,                     &
                            chunks(c)%field%celldx,                    &
                            chunks(c)%field%celldy,                    &
                            density0,                  &
                            pressure,                  &
                            viscosity,                 &
                            xvel0,                     &
                            yvel0                      )
      ENDIF

    ENDIF

  ENDDO

END SUBROUTINE viscosity_driver

END MODULE viscosity_module
