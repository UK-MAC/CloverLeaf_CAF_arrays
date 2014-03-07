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

!>  @brief Communication Utilities
!>  @author Wayne Gaudin
!>  @details Contains all utilities required to run CloverLeaf in a distributed
!>  environment, including initialisation, mesh decompostion, reductions and
!>  halo exchange using explicit buffers.
!>
!>  Note the halo exchange is currently coded as simply as possible and no 
!>  optimisations have been implemented, such as post receives before sends or packing
!>  buffers with multiple data fields. This is intentional so the effect of these
!>  optimisations can be measured on large systems, as and when they are added.
!>
!>  Even without these modifications CloverLeaf weak scales well on moderately sized
!>  systems of the order of 10K cores.

MODULE clover_module

  USE data_module
  USE definitions_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE clover_barrier

    sync all

END SUBROUTINE clover_barrier

SUBROUTINE clover_abort

    ERROR STOP

END SUBROUTINE clover_abort

SUBROUTINE clover_finalize

  INTEGER :: err

  CLOSE(g_out)
  CALL FLUSH(0)
  CALL FLUSH(6)
  CALL FLUSH(g_out)

END SUBROUTINE clover_finalize

SUBROUTINE clover_init_comms

  IMPLICIT NONE

  parallel%parallel=.TRUE.

  parallel%image = this_image()
  parallel%max_image = num_images()

  parallel%task = parallel%image - 1

  IF(parallel%image.EQ.1) THEN
    parallel%boss=.TRUE.
  ENDIF

  parallel%boss_task=0
  parallel%max_task=parallel%max_image

END SUBROUTINE clover_init_comms

SUBROUTINE clover_get_num_chunks(count)

  IMPLICIT NONE

  INTEGER :: count

! Should be changed so there can be more than one chunk per mpi task

  count=parallel%max_task

END SUBROUTINE clover_get_num_chunks

SUBROUTINE clover_decompose(x_cells,y_cells,left,right,bottom,top)

  ! This decomposes the mesh into a number of chunks.
  ! The number of chunks may be a multiple of the number of mpi tasks
  ! Doesn't always return the best split if there are few factors
  ! All factors need to be stored and the best picked. But its ok for now

  IMPLICIT NONE

  INTEGER :: x_cells,y_cells,left(:),right(:),top(:),bottom(:)
  INTEGER :: c,delta_x,delta_y

  REAL(KIND=8) :: mesh_ratio,factor_x,factor_y
  INTEGER  :: chunk_x,chunk_y,mod_x,mod_y,split_found

  INTEGER  :: cx,cy,chunk,add_x,add_y,add_x_prev,add_y_prev
#ifdef LOCAL_SYNC
  INTEGER :: numNeighbours,n
#endif

  ! 2D Decomposition of the mesh

  mesh_ratio=real(x_cells)/real(y_cells)

  chunk_x=number_of_chunks
  chunk_y=1

  split_found=0 ! Used to detect 1D decomposition
  DO c=1,number_of_chunks
    IF (MOD(number_of_chunks,c).EQ.0) THEN
      factor_x=number_of_chunks/real(c)
      factor_y=c
      !Compare the factor ratio with the mesh ratio
      IF(factor_x/factor_y.LE.mesh_ratio) THEN
        chunk_y=c
        chunk_x=number_of_chunks/c
        split_found=1
        EXIT
      ENDIF
    ENDIF
  ENDDO

  IF(split_found.EQ.0.OR.chunk_y.EQ.number_of_chunks) THEN ! Prime number or 1D decomp detected
    IF(mesh_ratio.GE.1.0) THEN
      chunk_x=number_of_chunks
      chunk_y=1
    ELSE
      chunk_x=1
      chunk_y=number_of_chunks
    ENDIF
  ENDIF

  delta_x=x_cells/chunk_x
  delta_y=y_cells/chunk_y
  mod_x=MOD(x_cells,chunk_x)
  mod_y=MOD(y_cells,chunk_y)

  ! Set up chunk mesh ranges and chunk connectivity

    add_x_prev=0
    add_y_prev=0
    chunk=1
    DO cy=1,chunk_y
        DO cx=1,chunk_x
            add_x=0
            add_y=0
            IF(cx.LE.mod_x)add_x=1
            IF(cy.LE.mod_y)add_y=1

            IF (chunk .EQ. parallel%task+1) THEN

                left(1)   = (cx-1)*delta_x+1+add_x_prev
                right(1)  = left(1)+delta_x-1+add_x
                bottom(1) = (cy-1)*delta_y+1+add_y_prev
                top(1)    = bottom(1)+delta_y-1+add_y

                chunks(1)%chunk_neighbours(chunk_left)=chunk_x*(cy-1)+cx-1
                chunks(1)%chunk_neighbours(chunk_right)=chunk_x*(cy-1)+cx+1
                chunks(1)%chunk_neighbours(chunk_bottom)=chunk_x*(cy-2)+cx
                chunks(1)%chunk_neighbours(chunk_top)=chunk_x*(cy)+cx

                IF(cx.EQ.1)       chunks(1)%chunk_neighbours(chunk_left)=external_face
                IF(cx.EQ.chunk_x) chunks(1)%chunk_neighbours(chunk_right)=external_face
                IF(cy.EQ.1)       chunks(1)%chunk_neighbours(chunk_bottom)=external_face
                IF(cy.EQ.chunk_y) chunks(1)%chunk_neighbours(chunk_top)=external_face

#ifdef LOCAL_SYNC
                numNeighbours=0
                IF (chunks(1)%chunk_neighbours(chunk_left).NE.external_face) THEN
                    numNeighbours = numNeighbours +1
                ENDIF
                IF (chunks(1)%chunk_neighbours(chunk_right).NE.external_face) THEN
                   numNeighbours = numNeighbours +1
                ENDIF
                IF (chunks(1)%chunk_neighbours(chunk_top).NE.external_face) THEN
                   numNeighbours = numNeighbours +1
                ENDIF
                IF (chunks(1)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
                   numNeighbours = numNeighbours +1
                ENDIF
                ALLOCATE(chunks(1)%imageNeighbours(numNeighbours))

                !caf:may need to update this when multiple chunks per image so that the image is recorded correctly 
                IF (numNeighbours > 0) THEN
                   n=1
                   IF (chunks(1)%chunk_neighbours(chunk_left).NE.external_face) THEN
                      chunks(1)%imageNeighbours(n) = chunks(1)%chunk_neighbours(chunk_left)
                      n=n+1
                   ENDIF
                   IF (chunks(1)%chunk_neighbours(chunk_right).NE.external_face) THEN
                      chunks(1)%imageNeighbours(n) = chunks(1)%chunk_neighbours(chunk_right)
                      n=n+1
                   ENDIF
                   IF (chunks(1)%chunk_neighbours(chunk_top).NE.external_face) THEN
                      chunks(1)%imageNeighbours(n) = chunks(1)%chunk_neighbours(chunk_top)
                      n=n+1
                   ENDIF
                   IF (chunks(1)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
                      chunks(1)%imageNeighbours(n) = chunks(1)%chunk_neighbours(chunk_bottom)
                      n=n+1
                   ENDIF
                ENDIF
#endif

            ENDIF

            IF(cx.LE.mod_x)add_x_prev=add_x_prev+1
            chunk=chunk+1
        ENDDO
        add_x_prev=0
        IF(cy.LE.mod_y)add_y_prev=add_y_prev+1
    ENDDO

  IF(parallel%boss)THEN
    WRITE(g_out,*)
    WRITE(g_out,*)"Mesh ratio of ",mesh_ratio
    WRITE(g_out,*)"Decomposing the mesh into ",chunk_x," by ",chunk_y," chunks"
    WRITE(g_out,*)
  ENDIF

END SUBROUTINE clover_decompose

SUBROUTINE clover_exchange(fields,depth)

    IMPLICIT NONE

    INTEGER      :: fields(:),depth
    INTEGER      :: chunk, left_neighbour, right_neighbour, top_neighbour, bottom_neighbour, xinc, yinc

    ! Assuming 1 patch per task, this will be changed
    ! Also, not packing all fields for each communication, doing one at a time

    ! the chunk number which the executing task is responsible for
    chunk = 1

    left_neighbour = chunks(chunk)%chunk_neighbours(chunk_left)
    right_neighbour = chunks(chunk)%chunk_neighbours(chunk_right)
    bottom_neighbour = chunks(chunk)%chunk_neighbours(chunk_bottom)
    top_neighbour = chunks(chunk)%chunk_neighbours(chunk_top)

    !caf: syncronise all processes here to ensure that all have finished the previous phase and there the halo exchange can begin
#ifdef LOCAL_SYNC
    sync images( chunks(chunk)%imageNeighbours )
#else
    sync all
#endif

    ! Assuming 1 patch per task, this will be changed
    ! Also, not packing all fields for each communication, doing one at a time

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        xinc = 0    
        yinc = 0

        IF(left_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the left of the current image
          density0(                                                 &
                   chunks(chunk)%field%x_max+xinc+1:chunks(chunk)%field%x_max+xinc+depth,      &
                   chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth )[left_neighbour] =    &
             density0(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                          chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF
        IF(right_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the right of the current image
          density0(                                              &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_min-1,              &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[right_neighbour] =   &
              density0(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,         &
                                           chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF

        ! caf: can be replaced with a sync with just the neighbours
#ifdef LOCAL_SYNC
        sync images( chunks(chunk)%imageNeighbours )
#else
        sync all
#endif

        IF(bottom_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          density0(                                             &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,    &
                  chunks(chunk)%field%y_max+yinc+1:chunks(chunk)%field%y_max+yinc+depth)[bottom_neighbour] = &
             density0(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                          chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
        ENDIF
        IF(top_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          density0(                                                  &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,         &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_min-1)[top_neighbour] =                &
             density0(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                          chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
        ENDIF
    ENDIF


    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        xinc = 0
        yinc = 0

        IF(left_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the left of the current image
          density1(                                                 &
                  chunks(chunk)%field%x_max+xinc+1:chunks(chunk)%field%x_max+xinc+depth,       &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[left_neighbour] =      &
             density1(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                          chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF
        IF(right_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the right of the current image
          density1(                                              &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_min-1,              &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[right_neighbour] =   &
             density1(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,          &
                                          chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF

        ! caf: can be replaced with a sync with just the neighbours
#ifdef LOCAL_SYNC
        sync images( chunks(chunk)%imageNeighbours )
#else
        sync all
#endif

        IF(bottom_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          density1(                                             &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,    &
                  chunks(chunk)%field%y_max+yinc+1:chunks(chunk)%field%y_max+yinc+depth)[bottom_neighbour] = &
             density1(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                          chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
        ENDIF
        IF(top_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          density1(                                                  &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,         &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_min-1)[top_neighbour] =                &
             density1(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                          chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
        ENDIF
    ENDIF


    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        xinc = 0
        yinc = 0

        IF(left_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the left of the current image
          energy0(                                                 &
                  chunks(chunk)%field%x_max+xinc+1:chunks(chunk)%field%x_max+xinc+depth,      &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[left_neighbour] =     &
             energy0(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                         chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF
        IF(right_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the right of the current image
          energy0(                                              &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_min-1,             &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[right_neighbour] =  &
             energy0(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,          &
                                         chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF

        ! caf: can be replaced with a sync with just the neighbours
#ifdef LOCAL_SYNC
        sync images( chunks(chunk)%imageNeighbours )
#else
        sync all
#endif

        IF(bottom_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          energy0(                                              &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,    &
                  chunks(chunk)%field%y_max+yinc+1:chunks(chunk)%field%y_max+yinc+depth)[bottom_neighbour] = &
             energy0(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,   &
                                         chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
        ENDIF
        IF(top_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          energy0(                                                  &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,        &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_min-1)[top_neighbour] =               &
             energy0(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                         chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
        ENDIF
    ENDIF


    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        xinc = 0
        yinc = 0

        IF(left_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the left of the current image
          energy1(                                                 &
                  chunks(chunk)%field%x_max+xinc+1:chunks(chunk)%field%x_max+xinc+depth,      &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[left_neighbour] =     &
             energy1(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                         chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF
        IF(right_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the right of the current image
          energy1(                                              &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_min-1,             &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[right_neighbour] =  &
             energy1(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,          &
                                         chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF

        ! caf: can be replaced with a sync with just the neighbours
#ifdef LOCAL_SYNC
        sync images( chunks(chunk)%imageNeighbours )
#else
        sync all
#endif

        IF(bottom_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          energy1(                                              &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,    &
                  chunks(chunk)%field%y_max+yinc+1:chunks(chunk)%field%y_max+yinc+depth)[bottom_neighbour] = &
             energy1(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,   &
                                         chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
        ENDIF
        IF(top_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          energy1(                                                  &
                   chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,       &
                   chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_min-1)[top_neighbour] =              &
             energy1(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                         chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
        ENDIF
    ENDIF


    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        xinc = 0
        yinc = 0

        IF(left_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the left of the current image
          pressure(                                                 &
                  chunks(chunk)%field%x_max+xinc+1:chunks(chunk)%field%x_max+xinc+depth,       &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[left_neighbour] =      &
             pressure(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                          chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF
        IF(right_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the right of the current image
          pressure(                                              &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_min-1,              &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[right_neighbour] =   &
             pressure(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,          &
                                          chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF

        ! caf: can be replaced with a sync with just the neighbours
#ifdef LOCAL_SYNC
        sync images( chunks(chunk)%imageNeighbours )
#else
        sync all
#endif

        IF(bottom_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          pressure(                                             &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,    &
                  chunks(chunk)%field%y_max+yinc+1:chunks(chunk)%field%y_max+yinc+depth)[bottom_neighbour] = &
             pressure(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                          chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
        ENDIF
        IF(top_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          pressure(                                                  &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,         &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_min-1)[top_neighbour] =                &
             pressure(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                          chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
        ENDIF
    ENDIF


    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        xinc = 0
        yinc = 0

        IF(left_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the left of the current image
          viscosity(                                                 &
                  chunks(chunk)%field%x_max+xinc+1:chunks(chunk)%field%x_max+xinc+depth,        &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[left_neighbour] =       &
             viscosity(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                           chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF
        IF(right_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the right of the current image
          viscosity(                                              &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_min-1,               &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[right_neighbour] =    &
             viscosity(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,          &
                                           chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF

        ! caf: can be replaced with a sync with just the neighbours
#ifdef LOCAL_SYNC
        sync images( chunks(chunk)%imageNeighbours )
#else
        sync all
#endif

        IF(bottom_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          viscosity(                                             &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,     &
                  chunks(chunk)%field%y_max+yinc+1:chunks(chunk)%field%y_max+yinc+depth)[bottom_neighbour] =  &
             viscosity(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                           chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
        ENDIF
        IF(top_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          viscosity(                                                   &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,           &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_min-1)[top_neighbour] =                  &
             viscosity(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                           chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
        ENDIF
    ENDIF


    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        xinc = 0
        yinc = 0

        IF(left_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the left of the current image
          soundspeed(                                                 &
                  chunks(chunk)%field%x_max+xinc+1:chunks(chunk)%field%x_max+xinc+depth,         &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[left_neighbour]=        &
             soundspeed(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                            chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF
        IF(right_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the right of the current image
          soundspeed(                                              &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_min-1,                &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[right_neighbour] =     &
             soundspeed(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,          &
                                            chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF

        ! caf: can be replaced with a sync with just the neighbours
#ifdef LOCAL_SYNC
        sync images( chunks(chunk)%imageNeighbours )
#else
        sync all
#endif

        IF(bottom_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          soundspeed(                                             &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,      &
                  chunks(chunk)%field%y_max+yinc+1:chunks(chunk)%field%y_max+yinc+depth)[bottom_neighbour] =   &
             soundspeed(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                            chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
        ENDIF
        IF(top_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          soundspeed(                                                   &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,            &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_min-1)[top_neighbour] =                   &
             soundspeed(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                            chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
        ENDIF
    ENDIF


    IF(fields(FIELD_XVEL0).EQ.1) THEN
        xinc = 1
        yinc = 1

        IF(left_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the left of the current image
          xvel0(                                                 &
                  chunks(chunk)%field%x_max+xinc+1:chunks(chunk)%field%x_max+xinc+depth,    &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[left_neighbour] =   &
             xvel0(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                       chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF
        IF(right_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the right of the current image
          xvel0(                                               &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_min-1,            &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[right_neighbour] = &
             xvel0(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,           &
                                       chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF

        ! caf: can be replaced with a sync with just the neighbours
#ifdef LOCAL_SYNC
        sync images( chunks(chunk)%imageNeighbours )
#else
        sync all
#endif

        IF(bottom_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          xvel0(                                                &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,    &
                  chunks(chunk)%field%y_max+yinc+1:chunks(chunk)%field%y_max+yinc+depth)[bottom_neighbour] = &
             xvel0(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,     &
                                       chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
        ENDIF
        IF(top_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          xvel0(                                                  &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,      &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_min-1)[top_neighbour] =             &
             xvel0(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                       chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
        ENDIF
    ENDIF


    IF(fields(FIELD_XVEL1).EQ.1) THEN
        xinc = 1
        yinc = 1

        IF(left_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the left of the current image
          xvel1(                                                 &
                  chunks(chunk)%field%x_max+xinc+1:chunks(chunk)%field%x_max+xinc+depth,    &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[left_neighbour] =   &
             xvel1(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                       chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF
        IF(right_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the right of the current image
          xvel1(                                               &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_min-1,            &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[right_neighbour] = &
             xvel1(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,           &
                                       chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF

        ! caf: can be replaced with a sync with just the neighbours
#ifdef LOCAL_SYNC
        sync images( chunks(chunk)%imageNeighbours )
#else
        sync all
#endif

        IF(bottom_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          xvel1(                                                &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,    &
                  chunks(chunk)%field%y_max+yinc+1:chunks(chunk)%field%y_max+yinc+depth)[bottom_neighbour] = &
             xvel1(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,     &
                                       chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
        ENDIF
        IF(top_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          xvel1(                                                  &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,      &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_min-1)[top_neighbour] =             &
             xvel1(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                       chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
        ENDIF
    ENDIF


    IF(fields(FIELD_YVEL0).EQ.1) THEN
        xinc = 1
        yinc = 1

        IF(left_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the left of the current image
          yvel0(                                                 &
                  chunks(chunk)%field%x_max+xinc+1:chunks(chunk)%field%x_max+xinc+depth,    &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[left_neighbour] =   &
             yvel0(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                       chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF
        IF(right_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the right of the current image
          yvel0(                                               &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_min-1,            &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[right_neighbour] = &
             yvel0(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,           &
                                       chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF

        ! caf: can be replaced with a sync with just the neighbours
#ifdef LOCAL_SYNC
        sync images( chunks(chunk)%imageNeighbours )
#else
        sync all
#endif

        IF(bottom_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          yvel0(                                                &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,    &
                  chunks(chunk)%field%y_max+yinc+1:chunks(chunk)%field%y_max+yinc+depth)[bottom_neighbour] = &
             yvel0(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,     &
                                       chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
        ENDIF
        IF(top_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          yvel0(                                                  &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,      &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_min-1)[top_neighbour] =             &
             yvel0(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                       chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
        ENDIF
    ENDIF


    IF(fields(FIELD_YVEL1).EQ.1) THEN
        xinc = 1
        yinc = 1

        IF(left_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the left of the current image
          yvel1(                                                 &
                  chunks(chunk)%field%x_max+xinc+1:chunks(chunk)%field%x_max+xinc+depth,    &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[left_neighbour] =   &
             yvel1(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                       chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF
        IF(right_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the right of the current image
          yvel1(                                               &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_min-1,            &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[right_neighbour] = &
             yvel1(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,           &
                                       chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF

        ! caf: can be replaced with a sync with just the neighbours
#ifdef LOCAL_SYNC
        sync images( chunks(chunk)%imageNeighbours )
#else
        sync all
#endif

        IF(bottom_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          yvel1(                                                &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,    &
                  chunks(chunk)%field%y_max+yinc+1:chunks(chunk)%field%y_max+yinc+depth)[bottom_neighbour] = &
             yvel1(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,     &
                                       chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
        ENDIF
        IF(top_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          yvel1(                                                  &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,      &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_min-1)[top_neighbour] =             &
             yvel1(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                       chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
        ENDIF
    ENDIF


    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        xinc = 1
        yinc = 0

        IF(left_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the left of the current image
          vol_flux_x(                                                 &
                  chunks(chunk)%field%x_max+xinc+1:chunks(chunk)%field%x_max+xinc+depth,         &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[left_neighbour] =        &
             vol_flux_x(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                            chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF
        IF(right_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the right of the current image
          vol_flux_x(                                              &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_min-1,                &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[right_neighbour] =     &
             vol_flux_x(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,          &
                                            chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF

        ! caf: can be replaced with a sync with just the neighbours
#ifdef LOCAL_SYNC
        sync images( chunks(chunk)%imageNeighbours )
#else
        sync all
#endif

        IF(bottom_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          vol_flux_x(                                             &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,      &
                  chunks(chunk)%field%y_max+yinc+1:chunks(chunk)%field%y_max+yinc+depth)[bottom_neighbour] =   &
             vol_flux_x(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                            chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
        ENDIF
        IF(top_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          vol_flux_x(                                                  &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,           &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_min-1)[top_neighbour] =                  &
             vol_flux_x(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                            chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
        ENDIF
    ENDIF


    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        xinc = 0
        yinc = 1

        IF(left_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the left of the current image
          vol_flux_y(                                                 &
                  chunks(chunk)%field%x_max+xinc+1:chunks(chunk)%field%x_max+xinc+depth,         &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[left_neighbour] =        &
             vol_flux_y(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                            chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF
        IF(right_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the right of the current image
          vol_flux_y(                                              &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_min-1,                &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[right_neighbour] =     &
             vol_flux_y(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,          &
                                            chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF

      ! caf: can be replaced with a sync with just the neighbours
#ifdef LOCAL_SYNC
        sync images( chunks(chunk)%imageNeighbours )
#else
        sync all
#endif

        IF(bottom_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          vol_flux_y(                                             &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,      &
                  chunks(chunk)%field%y_max+yinc+1:chunks(chunk)%field%y_max+yinc+depth)[bottom_neighbour] =   &
             vol_flux_y(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                            chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
        ENDIF
        IF(top_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          vol_flux_y(                                                  &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,           &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_min-1)[top_neighbour] =                  &
             vol_flux_y(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                            chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
        ENDIF
    ENDIF


    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        xinc = 1
        yinc = 0

        IF(left_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the left of the current image
          mass_flux_x(                                                 &
                  chunks(chunk)%field%x_max+xinc+1:chunks(chunk)%field%x_max+xinc+depth,          &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[left_neighbour] =         &
             mass_flux_x(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                             chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF
        IF(right_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the right of the current image
          mass_flux_x(                                              &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_min-1,                 &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[right_neighbour] =      &
             mass_flux_x(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,          &
                                             chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF

        ! caf: can be replaced with a sync with just the neighbours
#ifdef LOCAL_SYNC
        sync images( chunks(chunk)%imageNeighbours )
#else
        sync all
#endif

        IF(bottom_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          mass_flux_x(                                             &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,       &
                  chunks(chunk)%field%y_max+yinc+1:chunks(chunk)%field%y_max+yinc+depth)[bottom_neighbour] =    &
             mass_flux_x(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                             chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
        ENDIF
        IF(top_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          mass_flux_x(                                                  &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,            &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_min-1)[top_neighbour] =                   &
             mass_flux_x(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                             chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
        ENDIF
    ENDIF


    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        xinc = 0
        yinc = 1

        IF(left_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the left of the current image
          mass_flux_y(                                                 &
                  chunks(chunk)%field%x_max+xinc+1:chunks(chunk)%field%x_max+xinc+depth,          &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[left_neighbour] =         &
             mass_flux_y(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                             chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF
        IF(right_neighbour.NE.external_face) THEN
          !caf: one sided put to the image on the right of the current image
          mass_flux_y(                                              &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_min-1,                 &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)[right_neighbour] =      &
             mass_flux_y(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,          &
                                             chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
        ENDIF

        ! caf: can be replaced with a sync with just the neighbours
#ifdef LOCAL_SYNC
        sync images( chunks(chunk)%imageNeighbours )
#else
        sync all
#endif

        IF(bottom_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          mass_flux_y(                                             &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,       &
                  chunks(chunk)%field%y_max+yinc+1:chunks(chunk)%field%y_max+yinc+depth)[bottom_neighbour] =    &
             mass_flux_y(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                             chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
        ENDIF
        IF(top_neighbour.NE.external_face) THEN
          !caf: one sided put to the image under the current image
          mass_flux_y(                                                   &
                  chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,             &
                  chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_min-1)[top_neighbour] =                    &
             mass_flux_y(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                             chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
        ENDIF
    ENDIF

  ! caf: syncronise all processes here to ensure that they have all finished the halo exchange before they start the next phase 
#ifdef LOCAL_SYNC
    sync images( chunks(chunk)%imageNeighbours )
#else
    sync all
#endif

END SUBROUTINE clover_exchange

SUBROUTINE clover_sum(value)

  ! Only sums to the master

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: total

  total=value

  CALL CO_SUM(value, total, result_image=1)

  value=total

END SUBROUTINE clover_sum

SUBROUTINE clover_min(value)

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: minimum

  minimum=value

  CALL CO_MIN(value, minimum)

  value=minimum

END SUBROUTINE clover_min

SUBROUTINE clover_max(value)

    IMPLICIT NONE

    REAL(KIND=8) :: value

    REAL(KIND=8) :: maximum

    INTEGER :: err

    maximum=value

    CALL CO_MAX(value, maximum)

    value=maximum

END SUBROUTINE clover_max

SUBROUTINE clover_allgather(value,values)

    IMPLICIT NONE

    REAL(KIND=8) :: value

    REAL(KIND=8) :: values(parallel%max_task)

    INTEGER :: err

    values(1)=value ! Just to ensure it will work in serial

    totals(parallel%image)[1] = value

END SUBROUTINE clover_allgather

SUBROUTINE clover_check_error(error)

  IMPLICIT NONE

  INTEGER :: error

  INTEGER :: maximum

  maximum=error

  CALL CO_MAX(error, maximum)

  error=maximum

END SUBROUTINE clover_check_error


END MODULE clover_module
!=======
!  IMPLICIT NONE
!
!  INTEGER      :: fields(:),depth,location_of_tasks_chunk
!
!  ! Assuming 1 patch per task, this will be changed
!  ! Also, not packing all fields for each communication, doing one at a time
!
!    location_of_tasks_chunk = 1
!
!  IF(fields(FIELD_DENSITY0).EQ.1) THEN
!    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
!                                 chunks(location_of_tasks_chunk)%field%density0,      &
!                                 chunks(location_of_tasks_chunk)%left_snd_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%left_rcv_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%right_snd_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%right_rcv_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%bottom_snd_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%bottom_rcv_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%top_snd_buffer,                       &
!                                 chunks(location_of_tasks_chunk)%top_rcv_buffer,                       &
!                                 depth,CELL_DATA)
!  ENDIF
!
!  IF(fields(FIELD_DENSITY1).EQ.1) THEN
!    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
!                                 chunks(location_of_tasks_chunk)%field%density1,      &
!                                 chunks(location_of_tasks_chunk)%left_snd_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%left_rcv_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%right_snd_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%right_rcv_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%bottom_snd_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%bottom_rcv_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%top_snd_buffer,                       &
!                                 chunks(location_of_tasks_chunk)%top_rcv_buffer,                       &
!                                 depth,CELL_DATA)
!  ENDIF
!
!  IF(fields(FIELD_ENERGY0).EQ.1) THEN
!    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
!                                 chunks(location_of_tasks_chunk)%field%energy0,       &
!                                 chunks(location_of_tasks_chunk)%left_snd_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%left_rcv_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%right_snd_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%right_rcv_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%bottom_snd_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%bottom_rcv_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%top_snd_buffer,                       &
!                                 chunks(location_of_tasks_chunk)%top_rcv_buffer,                       &
!                                 depth,CELL_DATA)
!  ENDIF
!
!  IF(fields(FIELD_ENERGY1).EQ.1) THEN
!    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
!                                 chunks(location_of_tasks_chunk)%field%energy1,       &
!                                 chunks(location_of_tasks_chunk)%left_snd_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%left_rcv_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%right_snd_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%right_rcv_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%bottom_snd_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%bottom_rcv_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%top_snd_buffer,                       &
!                                 chunks(location_of_tasks_chunk)%top_rcv_buffer,                       &
!                                 depth,CELL_DATA)
!  ENDIF
!
!  IF(fields(FIELD_PRESSURE).EQ.1) THEN
!    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
!                                 chunks(location_of_tasks_chunk)%field%pressure,      &
!                                 chunks(location_of_tasks_chunk)%left_snd_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%left_rcv_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%right_snd_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%right_rcv_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%bottom_snd_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%bottom_rcv_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%top_snd_buffer,                       &
!                                 chunks(location_of_tasks_chunk)%top_rcv_buffer,                       &
!                                 depth,CELL_DATA)
!  ENDIF
!
!  IF(fields(FIELD_VISCOSITY).EQ.1) THEN
!    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
!                                 chunks(location_of_tasks_chunk)%field%viscosity,     &
!                                 chunks(location_of_tasks_chunk)%left_snd_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%left_rcv_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%right_snd_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%right_rcv_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%bottom_snd_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%bottom_rcv_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%top_snd_buffer,                       &
!                                 chunks(location_of_tasks_chunk)%top_rcv_buffer,                       &
!                                 depth,CELL_DATA)
!  ENDIF
!
!  IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
!    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
!                                 chunks(location_of_tasks_chunk)%field%soundspeed,    &
!                                 chunks(location_of_tasks_chunk)%left_snd_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%left_rcv_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%right_snd_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%right_rcv_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%bottom_snd_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%bottom_rcv_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%top_snd_buffer,                       &
!                                 chunks(location_of_tasks_chunk)%top_rcv_buffer,                       &
!                                 depth,CELL_DATA)
!  ENDIF
!
!  IF(fields(FIELD_XVEL0).EQ.1) THEN
!    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
!                                 chunks(location_of_tasks_chunk)%field%xvel0,         &
!                                 chunks(location_of_tasks_chunk)%left_snd_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%left_rcv_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%right_snd_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%right_rcv_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%bottom_snd_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%bottom_rcv_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%top_snd_buffer,                       &
!                                 chunks(location_of_tasks_chunk)%top_rcv_buffer,                       &
!                                 depth,VERTEX_DATA)
!  ENDIF
!
!  IF(fields(FIELD_XVEL1).EQ.1) THEN
!    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
!                                 chunks(location_of_tasks_chunk)%field%xvel1,         &
!                                 chunks(location_of_tasks_chunk)%left_snd_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%left_rcv_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%right_snd_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%right_rcv_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%bottom_snd_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%bottom_rcv_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%top_snd_buffer,                       &
!                                 chunks(location_of_tasks_chunk)%top_rcv_buffer,                       &
!                                 depth,VERTEX_DATA)
!  ENDIF
!
!  IF(fields(FIELD_YVEL0).EQ.1) THEN
!    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
!                                 chunks(location_of_tasks_chunk)%field%yvel0,         &
!                                 chunks(location_of_tasks_chunk)%left_snd_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%left_rcv_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%right_snd_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%right_rcv_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%bottom_snd_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%bottom_rcv_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%top_snd_buffer,                       &
!                                 chunks(location_of_tasks_chunk)%top_rcv_buffer,                       &
!                                 depth,VERTEX_DATA)
!  ENDIF
!
!  IF(fields(FIELD_YVEL1).EQ.1) THEN
!    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
!                                 chunks(location_of_tasks_chunk)%field%yvel1,         &
!                                 chunks(location_of_tasks_chunk)%left_snd_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%left_rcv_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%right_snd_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%right_rcv_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%bottom_snd_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%bottom_rcv_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%top_snd_buffer,                       &
!                                 chunks(location_of_tasks_chunk)%top_rcv_buffer,                       &
!                                 depth,VERTEX_DATA)
!  ENDIF
!
!  IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
!    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
!                                 chunks(location_of_tasks_chunk)%field%vol_flux_x,    &
!                                 chunks(location_of_tasks_chunk)%left_snd_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%left_rcv_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%right_snd_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%right_rcv_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%bottom_snd_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%bottom_rcv_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%top_snd_buffer,                       &
!                                 chunks(location_of_tasks_chunk)%top_rcv_buffer,                       &
!                                 depth,X_FACE_DATA)
!  ENDIF
!
!  IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
!    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
!                                 chunks(location_of_tasks_chunk)%field%vol_flux_y,    &
!                                 chunks(location_of_tasks_chunk)%left_snd_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%left_rcv_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%right_snd_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%right_rcv_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%bottom_snd_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%bottom_rcv_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%top_snd_buffer,                       &
!                                 chunks(location_of_tasks_chunk)%top_rcv_buffer,                       &
!                                 depth,Y_FACE_DATA)
!  ENDIF
!
!  IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
!    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
!                                 chunks(location_of_tasks_chunk)%field%mass_flux_x,   &
!                                 chunks(location_of_tasks_chunk)%left_snd_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%left_rcv_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%right_snd_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%right_rcv_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%bottom_snd_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%bottom_rcv_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%top_snd_buffer,                       &
!                                 chunks(location_of_tasks_chunk)%top_rcv_buffer,                       &
!                                 depth,X_FACE_DATA)
!  ENDIF
!
!  IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
!    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
!                                 chunks(location_of_tasks_chunk)%field%mass_flux_y,   &
!                                 chunks(location_of_tasks_chunk)%left_snd_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%left_rcv_buffer,                      &
!                                 chunks(location_of_tasks_chunk)%right_snd_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%right_rcv_buffer,                     &
!                                 chunks(location_of_tasks_chunk)%bottom_snd_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%bottom_rcv_buffer,                    &
!                                 chunks(location_of_tasks_chunk)%top_snd_buffer,                       &
!                                 chunks(location_of_tasks_chunk)%top_rcv_buffer,                       &
!                                 depth,Y_FACE_DATA)
!  ENDIF
!>>>>>>> 05155921477081faa03dc409886647d4ecacdc34
!=======
!    ! Send/receive the data
!    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
!      tag=4*(chunk)+1 ! 4 because we have 4 faces, 1 because it is leaving the left face
!      receiver=chunks(chunk)%chunk_neighbours(chunk_left)-1
!      CALL MPI_ISEND(left_snd_buffer,size,MPI_DOUBLE_PRECISION,receiver,tag &
!                    ,MPI_COMM_WORLD,request(message_count+1),err)
!
!      tag=4*(chunk)+2 ! 4 because we have 4 faces, 1 because it is coming from the right face of the left neighbour
!      sender=chunks(chunk)%chunk_neighbours(chunk_left)-1
!      CALL MPI_IRECV(left_rcv_buffer,size,MPI_DOUBLE_PRECISION,sender,tag &
!                    ,MPI_COMM_WORLD,request(message_count+2),err)
!      message_count=message_count+2
!    ENDIF
!
!    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
!      tag=4*chunk+2 ! 4 because we have 4 faces, 2 because it is leaving the right face
!      receiver=chunks(chunk)%chunk_neighbours(chunk_right)-1
!      CALL MPI_ISEND(right_snd_buffer,size,MPI_DOUBLE_PRECISION,receiver,tag &
!                    ,MPI_COMM_WORLD,request(message_count+1),err)
!
!      tag=4*(chunk)+1 ! 4 because we have 4 faces, 1 because it is coming from the left face of the right neighbour
!      sender=chunks(chunk)%chunk_neighbours(chunk_right)-1
!      CALL MPI_IRECV(right_rcv_buffer,size,MPI_DOUBLE_PRECISION,sender,tag, &
!                     MPI_COMM_WORLD,request(message_count+2),err)
!      message_count=message_count+2
!>>>>>>> 05155921477081faa03dc409886647d4ecacdc34
!=======
!    ! Send/receive the data
!    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
!      tag=4*(chunk)+3 ! 4 because we have 4 faces, 3 because it is leaving the bottom face
!      receiver=chunks(chunk)%chunk_neighbours(chunk_bottom)-1
!      CALL MPI_ISEND(bottom_snd_buffer,size,MPI_DOUBLE_PRECISION,receiver,tag &
!                    ,MPI_COMM_WORLD,request(message_count+1),err)
!
!      tag=4*(chunk)+4 ! 4 because we have 4 faces, 1 because it is coming from the top face of the bottom neighbour
!      sender=chunks(chunk)%chunk_neighbours(chunk_bottom)-1
!      CALL MPI_IRECV(bottom_rcv_buffer,size,MPI_DOUBLE_PRECISION,sender,tag &
!                    ,MPI_COMM_WORLD,request(message_count+2),err)
!      message_count=message_count+2
!    ENDIF
!
!    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
!      tag=4*(chunk)+4 ! 4 because we have 4 faces, 4 because it is leaving the top face
!      receiver=chunks(chunk)%chunk_neighbours(chunk_top)-1
!      CALL MPI_ISEND(top_snd_buffer,size,MPI_DOUBLE_PRECISION,receiver,tag &
!                    ,MPI_COMM_WORLD,request(message_count+1),err)
!
!      tag=4*(chunk)+3 ! 4 because we have 4 faces, 4 because it is coming from the left face of the top neighbour
!      sender=chunks(chunk)%chunk_neighbours(chunk_top)-1
!      CALL MPI_IRECV(top_rcv_buffer,size,MPI_DOUBLE_PRECISION,sender,tag, &
!                     MPI_COMM_WORLD,request(message_count+2),err)
!      message_count=message_count+2
!>>>>>>> 05155921477081faa03dc409886647d4ecacdc34
