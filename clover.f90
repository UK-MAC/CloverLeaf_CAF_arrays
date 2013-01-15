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
  USE MPI

  IMPLICIT NONE

CONTAINS

SUBROUTINE clover_barrier

  INTEGER :: err

  CALL MPI_BARRIER(MPI_COMM_WORLD,err)

END SUBROUTINE clover_barrier

SUBROUTINE clover_abort

  INTEGER :: ierr,err

  CALL MPI_ABORT(MPI_COMM_WORLD,ierr,err)

END SUBROUTINE clover_abort

SUBROUTINE clover_finalize

  INTEGER :: err

  CLOSE(g_out)
  CALL FLUSH(0)
  CALL FLUSH(6)
  CALL FLUSH(g_out)
  CALL MPI_FINALIZE(err)

END SUBROUTINE clover_finalize

SUBROUTINE clover_init_comms

  IMPLICIT NONE

  INTEGER :: err,rank,size

  rank=0
  size=1

  CALL MPI_INIT(err) 

  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,err) 
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,err) 

  parallel%parallel=.TRUE.
  parallel%task=rank

  IF(rank.EQ.0) THEN
    parallel%boss=.TRUE.
  ENDIF

  parallel%boss_task=0
  parallel%max_task=size

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
      left(chunk)=(cx-1)*delta_x+1+add_x_prev
      right(chunk)=left(chunk)+delta_x-1+add_x
      bottom(chunk)=(cy-1)*delta_y+1+add_y_prev
      top(chunk)=bottom(chunk)+delta_y-1+add_y
      chunks(chunk)%chunk_neighbours(chunk_left)=chunk_x*(cy-1)+cx-1
      chunks(chunk)%chunk_neighbours(chunk_right)=chunk_x*(cy-1)+cx+1
      chunks(chunk)%chunk_neighbours(chunk_bottom)=chunk_x*(cy-2)+cx
      chunks(chunk)%chunk_neighbours(chunk_top)=chunk_x*(cy)+cx
      IF(cx.EQ.1)chunks(chunk)%chunk_neighbours(chunk_left)=external_face
      IF(cx.EQ.chunk_x)chunks(chunk)%chunk_neighbours(chunk_right)=external_face
      IF(cy.EQ.1)chunks(chunk)%chunk_neighbours(chunk_bottom)=external_face
      IF(cy.EQ.chunk_y)chunks(chunk)%chunk_neighbours(chunk_top)=external_face
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

SUBROUTINE clover_allocate_buffers(chunk)

  IMPLICIT NONE

  INTEGER      :: chunk
  
  ! Unallocated buffers for external boundaries caused issues on some systems so they are now
  !  all allocated
  IF(parallel%task.EQ.chunks(chunk)%task)THEN
    !IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%left_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%left_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%right_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%right_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%bottom_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%bottom_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%top_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%top_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
    !ENDIF
  ENDIF

END SUBROUTINE clover_allocate_buffers

SUBROUTINE clover_exchange(fields,depth)

  IMPLICIT NONE

  INTEGER      :: fields(:),depth
  INTEGER      :: chunk, left_neighbour, right_neighbour, top_neighbour, bottom_neighbour, xinc, yinc

  ! Assuming 1 patch per task, this will be changed
  ! Also, not packing all fields for each communication, doing one at a time

  ! the chunk number which the executing task is responsible for
  chunk = parallel%task+1

  left_neighbour = chunks(chunk)%chunk_neighbours(chunk_left)
  right_neighbour = chunks(chunk)%chunk_neighbours(chunk_right)
  bottom_neighbour = chunks(chunk)%chunk_neighbours(chunk_bottom)
  top_neighbour = chunks(chunk)%chunk_neighbours(chunk_top)

  !caf: syncronise all processes here to ensure that all have finished the previous phase and there the halo exchange can begin
  sync all

  ! Assuming 1 patch per task, this will be changed
  ! Also, not packing all fields for each communication, doing one at a time

  IF(fields(FIELD_DENSITY0).EQ.1) THEN
    xinc = 0    
    yinc = 0

    IF(left_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the left of the current image
      chunks(left_neighbour)[left_neighbour]%field%density0(                                                 &
               chunks(left_neighbour)%field%x_max+xinc+1:chunks(left_neighbour)%field%x_max+xinc+depth,      &
               chunks(left_neighbour)%field%y_min-depth:chunks(left_neighbour)%field%y_max+yinc+depth ) =    &
         chunks(chunk)%field%density0(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                      chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF
    IF(right_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the right of the current image
      chunks(right_neighbour)[right_neighbour]%field%density0(                                              &
              chunks(right_neighbour)%field%x_min-depth:chunks(right_neighbour)%field%x_min-1,              &
              chunks(right_neighbour)%field%y_min-depth:chunks(right_neighbour)%field%y_max+yinc+depth) =   &
          chunks(chunk)%field%density0(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,         &
                                       chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF

    ! caf: can be replaced with a sync with just the neighbours
    sync all

    IF(bottom_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(bottom_neighbour)[bottom_neighbour]%field%density0(                                             &
              chunks(bottom_neighbour)%field%x_min-depth:chunks(bottom_neighbour)%field%x_max+xinc+depth,    &
              chunks(bottom_neighbour)%field%y_max+yinc+1:chunks(bottom_neighbour)%field%y_max+yinc+depth) = &
         chunks(chunk)%field%density0(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                      chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
    ENDIF
    IF(top_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(top_neighbour)[top_neighbour]%field%density0(                                                  &
              chunks(top_neighbour)%field%x_min-depth:chunks(top_neighbour)%field%x_max+xinc+depth,         &
              chunks(top_neighbour)%field%y_min-depth:chunks(top_neighbour)%field%y_min-1) =                &
         chunks(chunk)%field%density0(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                      chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
    ENDIF
  ENDIF

  IF(fields(FIELD_DENSITY1).EQ.1) THEN
    xinc = 0
    yinc = 0

    IF(left_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the left of the current image
      chunks(left_neighbour)[left_neighbour]%field%density1(                                                 &
              chunks(left_neighbour)%field%x_max+xinc+1:chunks(left_neighbour)%field%x_max+xinc+depth,       &
              chunks(left_neighbour)%field%y_min-depth:chunks(left_neighbour)%field%y_max+yinc+depth) =      &
         chunks(chunk)%field%density1(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                      chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF
    IF(right_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the right of the current image
      chunks(right_neighbour)[right_neighbour]%field%density1(                                              &
              chunks(right_neighbour)%field%x_min-depth:chunks(right_neighbour)%field%x_min-1,              &
              chunks(right_neighbour)%field%y_min-depth:chunks(right_neighbour)%field%y_max+yinc+depth) =   &
         chunks(chunk)%field%density1(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,          &
                                      chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF

    ! caf: can be replaced with a sync with just the neighbours
    sync all

    IF(bottom_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(bottom_neighbour)[bottom_neighbour]%field%density1(                                             &
              chunks(bottom_neighbour)%field%x_min-depth:chunks(bottom_neighbour)%field%x_max+xinc+depth,    &
              chunks(bottom_neighbour)%field%y_max+yinc+1:chunks(bottom_neighbour)%field%y_max+yinc+depth) = &
         chunks(chunk)%field%density1(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                      chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
    ENDIF
    IF(top_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(top_neighbour)[top_neighbour]%field%density1(                                                  &
              chunks(top_neighbour)%field%x_min-depth:chunks(top_neighbour)%field%x_max+xinc+depth,         &
              chunks(top_neighbour)%field%y_min-depth:chunks(top_neighbour)%field%y_min-1) =                &
         chunks(chunk)%field%density1(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                      chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
    ENDIF
  ENDIF

  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    xinc = 0
    yinc = 0

    IF(left_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the left of the current image
      chunks(left_neighbour)[left_neighbour]%field%energy0(                                                 &
              chunks(left_neighbour)%field%x_max+xinc+1:chunks(left_neighbour)%field%x_max+xinc+depth,      &
              chunks(left_neighbour)%field%y_min-depth:chunks(left_neighbour)%field%y_max+yinc+depth) =     &
         chunks(chunk)%field%energy0(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                     chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF
    IF(right_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the right of the current image
      chunks(right_neighbour)[right_neighbour]%field%energy0(                                              &
              chunks(right_neighbour)%field%x_min-depth:chunks(right_neighbour)%field%x_min-1,             &
              chunks(right_neighbour)%field%y_min-depth:chunks(right_neighbour)%field%y_max+yinc+depth) =  &
         chunks(chunk)%field%energy0(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,          &
                                     chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF

    ! caf: can be replaced with a sync with just the neighbours
    sync all

    IF(bottom_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(bottom_neighbour)[bottom_neighbour]%field%energy0(                                              &
              chunks(bottom_neighbour)%field%x_min-depth:chunks(bottom_neighbour)%field%x_max+xinc+depth,    &
              chunks(bottom_neighbour)%field%y_max+yinc+1:chunks(bottom_neighbour)%field%y_max+yinc+depth) = &
         chunks(chunk)%field%energy0(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,   &
                                     chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
    ENDIF
    IF(top_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(top_neighbour)[top_neighbour]%field%energy0(                                                  &
              chunks(top_neighbour)%field%x_min-depth:chunks(top_neighbour)%field%x_max+xinc+depth,        &
              chunks(top_neighbour)%field%y_min-depth:chunks(top_neighbour)%field%y_min-1) =               &
         chunks(chunk)%field%energy0(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                     chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
    ENDIF
  ENDIF

  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    xinc = 0
    yinc = 0

    IF(left_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the left of the current image
      chunks(left_neighbour)[left_neighbour]%field%energy1(                                                 &
              chunks(left_neighbour)%field%x_max+xinc+1:chunks(left_neighbour)%field%x_max+xinc+depth,      &
              chunks(left_neighbour)%field%y_min-depth:chunks(left_neighbour)%field%y_max+yinc+depth) =     &
         chunks(chunk)%field%energy1(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                     chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF
    IF(right_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the right of the current image
      chunks(right_neighbour)[right_neighbour]%field%energy1(                                              &
              chunks(right_neighbour)%field%x_min-depth:chunks(right_neighbour)%field%x_min-1,             &
              chunks(right_neighbour)%field%y_min-depth:chunks(right_neighbour)%field%y_max+yinc+depth) =  &
         chunks(chunk)%field%energy1(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,          &
                                     chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF

    ! caf: can be replaced with a sync with just the neighbours
    sync all

    IF(bottom_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(bottom_neighbour)[bottom_neighbour]%field%energy1(                                              &
              chunks(bottom_neighbour)%field%x_min-depth:chunks(bottom_neighbour)%field%x_max+xinc+depth,    &
              chunks(bottom_neighbour)%field%y_max+yinc+1:chunks(bottom_neighbour)%field%y_max+yinc+depth) = &
         chunks(chunk)%field%energy1(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,   &
                                     chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
    ENDIF
    IF(top_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(top_neighbour)[top_neighbour]%field%energy1(                                                  &
               chunks(top_neighbour)%field%x_min-depth:chunks(top_neighbour)%field%x_max+xinc+depth,       &
               chunks(top_neighbour)%field%y_min-depth:chunks(top_neighbour)%field%y_min-1) =              &
         chunks(chunk)%field%energy1(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                     chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
    ENDIF
  ENDIF

  IF(fields(FIELD_PRESSURE).EQ.1) THEN
    xinc = 0
    yinc = 0

    IF(left_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the left of the current image
      chunks(left_neighbour)[left_neighbour]%field%pressure(                                                 &
              chunks(left_neighbour)%field%x_max+xinc+1:chunks(left_neighbour)%field%x_max+xinc+depth,       &
              chunks(left_neighbour)%field%y_min-depth:chunks(left_neighbour)%field%y_max+yinc+depth) =      &
         chunks(chunk)%field%pressure(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                      chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF
    IF(right_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the right of the current image
      chunks(right_neighbour)[right_neighbour]%field%pressure(                                              &
              chunks(right_neighbour)%field%x_min-depth:chunks(right_neighbour)%field%x_min-1,              &
              chunks(right_neighbour)%field%y_min-depth:chunks(right_neighbour)%field%y_max+yinc+depth) =   &
         chunks(chunk)%field%pressure(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,          &
                                      chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF

    ! caf: can be replaced with a sync with just the neighbours
    sync all

    IF(bottom_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(bottom_neighbour)[bottom_neighbour]%field%pressure(                                             &
              chunks(bottom_neighbour)%field%x_min-depth:chunks(bottom_neighbour)%field%x_max+xinc+depth,    &
              chunks(bottom_neighbour)%field%y_max+yinc+1:chunks(bottom_neighbour)%field%y_max+yinc+depth) = &
         chunks(chunk)%field%pressure(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                      chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
    ENDIF
    IF(top_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(top_neighbour)[top_neighbour]%field%pressure(                                                  &
              chunks(top_neighbour)%field%x_min-depth:chunks(top_neighbour)%field%x_max+xinc+depth,         &
              chunks(top_neighbour)%field%y_min-depth:chunks(top_neighbour)%field%y_min-1) =                &
         chunks(chunk)%field%pressure(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                      chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
    ENDIF
  ENDIF

  IF(fields(FIELD_VISCOSITY).EQ.1) THEN
    xinc = 0
    yinc = 0

    IF(left_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the left of the current image
      chunks(left_neighbour)[left_neighbour]%field%viscosity(                                                 &
              chunks(left_neighbour)%field%x_max+xinc+1:chunks(left_neighbour)%field%x_max+xinc+depth,        &
              chunks(left_neighbour)%field%y_min-depth:chunks(left_neighbour)%field%y_max+yinc+depth) =       &
         chunks(chunk)%field%viscosity(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                       chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF
    IF(right_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the right of the current image
      chunks(right_neighbour)[right_neighbour]%field%viscosity(                                              &
              chunks(right_neighbour)%field%x_min-depth:chunks(right_neighbour)%field%x_min-1,               &
              chunks(right_neighbour)%field%y_min-depth:chunks(right_neighbour)%field%y_max+yinc+depth) =    &
         chunks(chunk)%field%viscosity(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,          &
                                       chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF

    ! caf: can be replaced with a sync with just the neighbours
    sync all

    IF(bottom_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(bottom_neighbour)[bottom_neighbour]%field%viscosity(                                             &
              chunks(bottom_neighbour)%field%x_min-depth:chunks(bottom_neighbour)%field%x_max+xinc+depth,     &
              chunks(bottom_neighbour)%field%y_max+yinc+1:chunks(bottom_neighbour)%field%y_max+yinc+depth) =  &
         chunks(chunk)%field%viscosity(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                       chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
    ENDIF
    IF(top_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(top_neighbour)[top_neighbour]%field%viscosity(                                                   &
              chunks(top_neighbour)%field%x_min-depth:chunks(top_neighbour)%field%x_max+xinc+depth,           &
              chunks(top_neighbour)%field%y_min-depth:chunks(top_neighbour)%field%y_min-1) =                  &
         chunks(chunk)%field%viscosity(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                       chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
    ENDIF
  ENDIF

  IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
    xinc = 0
    yinc = 0

    IF(left_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the left of the current image
      chunks(left_neighbour)[left_neighbour]%field%soundspeed(                                                 &
              chunks(left_neighbour)%field%x_max+xinc+1:chunks(left_neighbour)%field%x_max+xinc+depth,         &
              chunks(left_neighbour)%field%y_min-depth:chunks(left_neighbour)%field%y_max+yinc+depth) =        &
         chunks(chunk)%field%soundspeed(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                        chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF
    IF(right_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the right of the current image
      chunks(right_neighbour)[right_neighbour]%field%soundspeed(                                              &
              chunks(right_neighbour)%field%x_min-depth:chunks(right_neighbour)%field%x_min-1,                &
              chunks(right_neighbour)%field%y_min-depth:chunks(right_neighbour)%field%y_max+yinc+depth) =     &
         chunks(chunk)%field%soundspeed(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,          &
                                        chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF

    ! caf: can be replaced with a sync with just the neighbours
    sync all

    IF(bottom_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(bottom_neighbour)[bottom_neighbour]%field%soundspeed(                                             &
              chunks(bottom_neighbour)%field%x_min-depth:chunks(bottom_neighbour)%field%x_max+xinc+depth,      &
              chunks(bottom_neighbour)%field%y_max+yinc+1:chunks(bottom_neighbour)%field%y_max+yinc+depth) =   &
         chunks(chunk)%field%soundspeed(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                        chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
    ENDIF
    IF(top_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(top_neighbour)[top_neighbour]%field%soundspeed(                                                   &
              chunks(top_neighbour)%field%x_min-depth:chunks(top_neighbour)%field%x_max+xinc+depth,            &
              chunks(top_neighbour)%field%y_min-depth:chunks(top_neighbour)%field%y_min-1) =                   &
         chunks(chunk)%field%soundspeed(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                        chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
    ENDIF
  ENDIF

  IF(fields(FIELD_XVEL0).EQ.1) THEN
    xinc = 1
    yinc = 1

    IF(left_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the left of the current image
      chunks(left_neighbour)[left_neighbour]%field%xvel0(                                                 &
              chunks(left_neighbour)%field%x_max+xinc+1:chunks(left_neighbour)%field%x_max+xinc+depth,    &
              chunks(left_neighbour)%field%y_min-depth:chunks(left_neighbour)%field%y_max+yinc+depth) =   &
         chunks(chunk)%field%xvel0(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                   chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF
    IF(right_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the right of the current image
      chunks(right_neighbour)[right_neighbour]%field%xvel0(                                               &
              chunks(right_neighbour)%field%x_min-depth:chunks(right_neighbour)%field%x_min-1,            &
              chunks(right_neighbour)%field%y_min-depth:chunks(right_neighbour)%field%y_max+yinc+depth) = &
         chunks(chunk)%field%xvel0(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,           &
                                   chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF

    ! caf: can be replaced with a sync with just the neighbours
    sync all

    IF(bottom_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(bottom_neighbour)[bottom_neighbour]%field%xvel0(                                                &
              chunks(bottom_neighbour)%field%x_min-depth:chunks(bottom_neighbour)%field%x_max+xinc+depth,    &
              chunks(bottom_neighbour)%field%y_max+yinc+1:chunks(bottom_neighbour)%field%y_max+yinc+depth) = &
         chunks(chunk)%field%xvel0(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,     &
                                   chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
    ENDIF
    IF(top_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(top_neighbour)[top_neighbour]%field%xvel0(                                                  &
              chunks(top_neighbour)%field%x_min-depth:chunks(top_neighbour)%field%x_max+xinc+depth,      &
              chunks(top_neighbour)%field%y_min-depth:chunks(top_neighbour)%field%y_min-1) =             &
         chunks(chunk)%field%xvel0(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                   chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
    ENDIF
  ENDIF

  IF(fields(FIELD_XVEL1).EQ.1) THEN
    xinc = 1
    yinc = 1

    IF(left_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the left of the current image
      chunks(left_neighbour)[left_neighbour]%field%xvel1(                                                 &
              chunks(left_neighbour)%field%x_max+xinc+1:chunks(left_neighbour)%field%x_max+xinc+depth,    &
              chunks(left_neighbour)%field%y_min-depth:chunks(left_neighbour)%field%y_max+yinc+depth) =   &
         chunks(chunk)%field%xvel1(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                   chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF
    IF(right_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the right of the current image
      chunks(right_neighbour)[right_neighbour]%field%xvel1(                                               &
              chunks(right_neighbour)%field%x_min-depth:chunks(right_neighbour)%field%x_min-1,            &
              chunks(right_neighbour)%field%y_min-depth:chunks(right_neighbour)%field%y_max+yinc+depth) = &
         chunks(chunk)%field%xvel1(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,           &
                                   chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF

    ! caf: can be replaced with a sync with just the neighbours
    sync all

    IF(bottom_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(bottom_neighbour)[bottom_neighbour]%field%xvel1(                                                &
              chunks(bottom_neighbour)%field%x_min-depth:chunks(bottom_neighbour)%field%x_max+xinc+depth,    &
              chunks(bottom_neighbour)%field%y_max+yinc+1:chunks(bottom_neighbour)%field%y_max+yinc+depth) = &
         chunks(chunk)%field%xvel1(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,     &
                                   chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
    ENDIF
    IF(top_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(top_neighbour)[top_neighbour]%field%xvel1(                                                  &
              chunks(top_neighbour)%field%x_min-depth:chunks(top_neighbour)%field%x_max+xinc+depth,      &
              chunks(top_neighbour)%field%y_min-depth:chunks(top_neighbour)%field%y_min-1) =             &
         chunks(chunk)%field%xvel1(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                   chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
    ENDIF
  ENDIF

  IF(fields(FIELD_YVEL0).EQ.1) THEN
    xinc = 1
    yinc = 1

    IF(left_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the left of the current image
      chunks(left_neighbour)[left_neighbour]%field%yvel0(                                                 &
              chunks(left_neighbour)%field%x_max+xinc+1:chunks(left_neighbour)%field%x_max+xinc+depth,    &
              chunks(left_neighbour)%field%y_min-depth:chunks(left_neighbour)%field%y_max+yinc+depth) =   &
         chunks(chunk)%field%yvel0(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                   chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF
    IF(right_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the right of the current image
      chunks(right_neighbour)[right_neighbour]%field%yvel0(                                               &
              chunks(right_neighbour)%field%x_min-depth:chunks(right_neighbour)%field%x_min-1,            &
              chunks(right_neighbour)%field%y_min-depth:chunks(right_neighbour)%field%y_max+yinc+depth) = &
         chunks(chunk)%field%yvel0(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,           &
                                   chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF

    ! caf: can be replaced with a sync with just the neighbours
    sync all

    IF(bottom_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(bottom_neighbour)[bottom_neighbour]%field%yvel0(                                                &
              chunks(bottom_neighbour)%field%x_min-depth:chunks(bottom_neighbour)%field%x_max+xinc+depth,    &
              chunks(bottom_neighbour)%field%y_max+yinc+1:chunks(bottom_neighbour)%field%y_max+yinc+depth) = &
         chunks(chunk)%field%yvel0(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,     &
                                   chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
    ENDIF
    IF(top_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(top_neighbour)[top_neighbour]%field%yvel0(                                                  &
              chunks(top_neighbour)%field%x_min-depth:chunks(top_neighbour)%field%x_max+xinc+depth,      &
              chunks(top_neighbour)%field%y_min-depth:chunks(top_neighbour)%field%y_min-1) =             &
         chunks(chunk)%field%yvel0(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                   chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
    ENDIF
  ENDIF

  IF(fields(FIELD_YVEL1).EQ.1) THEN
    xinc = 1
    yinc = 1

    IF(left_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the left of the current image
      chunks(left_neighbour)[left_neighbour]%field%yvel1(                                                 &
              chunks(left_neighbour)%field%x_max+xinc+1:chunks(left_neighbour)%field%x_max+xinc+depth,    &
              chunks(left_neighbour)%field%y_min-depth:chunks(left_neighbour)%field%y_max+yinc+depth) =   &
         chunks(chunk)%field%yvel1(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                   chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF
    IF(right_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the right of the current image
      chunks(right_neighbour)[right_neighbour]%field%yvel1(                                               &
              chunks(right_neighbour)%field%x_min-depth:chunks(right_neighbour)%field%x_min-1,            &
              chunks(right_neighbour)%field%y_min-depth:chunks(right_neighbour)%field%y_max+yinc+depth) = &
         chunks(chunk)%field%yvel1(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,           &
                                   chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF

    ! caf: can be replaced with a sync with just the neighbours
    sync all

    IF(bottom_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(bottom_neighbour)[bottom_neighbour]%field%yvel1(                                                &
              chunks(bottom_neighbour)%field%x_min-depth:chunks(bottom_neighbour)%field%x_max+xinc+depth,    &
              chunks(bottom_neighbour)%field%y_max+yinc+1:chunks(bottom_neighbour)%field%y_max+yinc+depth) = &
         chunks(chunk)%field%yvel1(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,     &
                                   chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
    ENDIF
    IF(top_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(top_neighbour)[top_neighbour]%field%yvel1(                                                  &
              chunks(top_neighbour)%field%x_min-depth:chunks(top_neighbour)%field%x_max+xinc+depth,      &
              chunks(top_neighbour)%field%y_min-depth:chunks(top_neighbour)%field%y_min-1) =             &
         chunks(chunk)%field%yvel1(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                   chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
    ENDIF
  ENDIF

  IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
    xinc = 1
    yinc = 0

    IF(left_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the left of the current image
      chunks(left_neighbour)[left_neighbour]%field%vol_flux_x(                                                 &
              chunks(left_neighbour)%field%x_max+xinc+1:chunks(left_neighbour)%field%x_max+xinc+depth,         &
              chunks(left_neighbour)%field%y_min-depth:chunks(left_neighbour)%field%y_max+yinc+depth) =        &
         chunks(chunk)%field%vol_flux_x(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                        chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF
    IF(right_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the right of the current image
      chunks(right_neighbour)[right_neighbour]%field%vol_flux_x(                                              &
              chunks(right_neighbour)%field%x_min-depth:chunks(right_neighbour)%field%x_min-1,                &
              chunks(right_neighbour)%field%y_min-depth:chunks(right_neighbour)%field%y_max+yinc+depth) =     &
         chunks(chunk)%field%vol_flux_x(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,          &
                                        chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF

    ! caf: can be replaced with a sync with just the neighbours
    sync all

    IF(bottom_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(bottom_neighbour)[bottom_neighbour]%field%vol_flux_x(                                             &
              chunks(bottom_neighbour)%field%x_min-depth:chunks(bottom_neighbour)%field%x_max+xinc+depth,      &
              chunks(bottom_neighbour)%field%y_max+yinc+1:chunks(bottom_neighbour)%field%y_max+yinc+depth) =   &
         chunks(chunk)%field%vol_flux_x(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                        chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
    ENDIF
    IF(top_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(top_neighbour)[top_neighbour]%field%vol_flux_x(                                                  &
              chunks(top_neighbour)%field%x_min-depth:chunks(top_neighbour)%field%x_max+xinc+depth,           &
              chunks(top_neighbour)%field%y_min-depth:chunks(top_neighbour)%field%y_min-1) =                  &
         chunks(chunk)%field%vol_flux_x(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                        chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
    ENDIF
  ENDIF

  IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
    xinc = 0
    yinc = 1

    IF(left_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the left of the current image
      chunks(left_neighbour)[left_neighbour]%field%vol_flux_y(                                                 &
              chunks(left_neighbour)%field%x_max+xinc+1:chunks(left_neighbour)%field%x_max+xinc+depth,         &
              chunks(left_neighbour)%field%y_min-depth:chunks(left_neighbour)%field%y_max+yinc+depth) =        &
         chunks(chunk)%field%vol_flux_y(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                        chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF
    IF(right_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the right of the current image
      chunks(right_neighbour)[right_neighbour]%field%vol_flux_y(                                              &
              chunks(right_neighbour)%field%x_min-depth:chunks(right_neighbour)%field%x_min-1,                &
              chunks(right_neighbour)%field%y_min-depth:chunks(right_neighbour)%field%y_max+yinc+depth) =     &
         chunks(chunk)%field%vol_flux_y(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,          &
                                        chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF

    ! caf: can be replaced with a sync with just the neighbours
    sync all

    IF(bottom_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(bottom_neighbour)[bottom_neighbour]%field%vol_flux_y(                                             &
              chunks(bottom_neighbour)%field%x_min-depth:chunks(bottom_neighbour)%field%x_max+xinc+depth,      &
              chunks(bottom_neighbour)%field%y_max+yinc+1:chunks(bottom_neighbour)%field%y_max+yinc+depth) =   &
         chunks(chunk)%field%vol_flux_y(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                        chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
    ENDIF
    IF(top_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(top_neighbour)[top_neighbour]%field%vol_flux_y(                                                  &
              chunks(top_neighbour)%field%x_min-depth:chunks(top_neighbour)%field%x_max+xinc+depth,           &
              chunks(top_neighbour)%field%y_min-depth:chunks(top_neighbour)%field%y_min-1) =                  &
         chunks(chunk)%field%vol_flux_y(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                        chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
    ENDIF
  ENDIF

  IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
    xinc = 1
    yinc = 0

    IF(left_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the left of the current image
      chunks(left_neighbour)[left_neighbour]%field%mass_flux_x(                                                 &
              chunks(left_neighbour)%field%x_max+xinc+1:chunks(left_neighbour)%field%x_max+xinc+depth,          &
              chunks(left_neighbour)%field%y_min-depth:chunks(left_neighbour)%field%y_max+yinc+depth) =         &
         chunks(chunk)%field%mass_flux_x(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                         chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF
    IF(right_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the right of the current image
      chunks(right_neighbour)[right_neighbour]%field%mass_flux_x(                                              &
              chunks(right_neighbour)%field%x_min-depth:chunks(right_neighbour)%field%x_min-1,                 &
              chunks(right_neighbour)%field%y_min-depth:chunks(right_neighbour)%field%y_max+yinc+depth) =      &
         chunks(chunk)%field%mass_flux_x(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,          &
                                         chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF

    ! caf: can be replaced with a sync with just the neighbours
    sync all

    IF(bottom_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(bottom_neighbour)[bottom_neighbour]%field%mass_flux_x(                                             &
              chunks(bottom_neighbour)%field%x_min-depth:chunks(bottom_neighbour)%field%x_max+xinc+depth,       &
              chunks(bottom_neighbour)%field%y_max+yinc+1:chunks(bottom_neighbour)%field%y_max+yinc+depth) =    &
         chunks(chunk)%field%mass_flux_x(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                         chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
    ENDIF
    IF(top_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(top_neighbour)[top_neighbour]%field%mass_flux_x(                                                  &
              chunks(top_neighbour)%field%x_min-depth:chunks(top_neighbour)%field%x_max+xinc+depth,            &
              chunks(top_neighbour)%field%y_min-depth:chunks(top_neighbour)%field%y_min-1) =                   &
         chunks(chunk)%field%mass_flux_x(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth, &
                                         chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
    ENDIF
  ENDIF

  IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
    xinc = 0
    yinc = 1

    IF(left_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the left of the current image
      chunks(left_neighbour)[left_neighbour]%field%mass_flux_y(                                                 &
              chunks(left_neighbour)%field%x_max+xinc+1:chunks(left_neighbour)%field%x_max+xinc+depth,          &
              chunks(left_neighbour)%field%y_min-depth:chunks(left_neighbour)%field%y_max+yinc+depth) =         &
         chunks(chunk)%field%mass_flux_y(chunks(chunk)%field%x_min+xinc:chunks(chunk)%field%x_min+xinc-1+depth, &
                                         chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF
    IF(right_neighbour.NE.external_face) THEN
      !caf: one sided put to the image on the right of the current image
      chunks(right_neighbour)[right_neighbour]%field%mass_flux_y(                                              &
              chunks(right_neighbour)%field%x_min-depth:chunks(right_neighbour)%field%x_min-1,                 &
              chunks(right_neighbour)%field%y_min-depth:chunks(right_neighbour)%field%y_max+yinc+depth) =      &
         chunks(chunk)%field%mass_flux_y(chunks(chunk)%field%x_max+1-depth:chunks(chunk)%field%x_max,          &
                                         chunks(chunk)%field%y_min-depth:chunks(chunk)%field%y_max+yinc+depth)
    ENDIF

    ! caf: can be replaced with a sync with just the neighbours
    sync all

    IF(bottom_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(bottom_neighbour)[bottom_neighbour]%field%mass_flux_y(                                             &
              chunks(bottom_neighbour)%field%x_min-depth:chunks(bottom_neighbour)%field%x_max+xinc+depth,       &
              chunks(bottom_neighbour)%field%y_max+yinc+1:chunks(bottom_neighbour)%field%y_max+yinc+depth) =    &
         chunks(chunk)%field%mass_flux_y(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                         chunks(chunk)%field%y_min+yinc:chunks(chunk)%field%y_min+yinc-1+depth)
    ENDIF
    IF(top_neighbour.NE.external_face) THEN
      !caf: one sided put to the image under the current image
      chunks(top_neighbour)[top_neighbour]%field%mass_flux_y(                                                   &
              chunks(top_neighbour)%field%x_min-depth:chunks(top_neighbour)%field%x_max+xinc+depth,             &
              chunks(top_neighbour)%field%y_min-depth:chunks(top_neighbour)%field%y_min-1) =                    &
         chunks(chunk)%field%mass_flux_y(chunks(chunk)%field%x_min-depth:chunks(chunk)%field%x_max+xinc+depth,  &
                                         chunks(chunk)%field%y_max+1-depth:chunks(chunk)%field%y_max)
    ENDIF
  ENDIF

  ! caf: syncronise all processes here to ensure that they have all finished the halo exchange before they start the next phase 
  sync all

END SUBROUTINE clover_exchange

SUBROUTINE clover_sum(value)

  ! Only sums to the master

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: total

  INTEGER :: err

  total=value

  CALL MPI_REDUCE(value,total,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,err)

  value=total

END SUBROUTINE clover_sum

SUBROUTINE clover_min(value)

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: minimum

  INTEGER :: err

  minimum=value

  CALL MPI_ALLREDUCE(value,minimum,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,err)

  value=minimum

END SUBROUTINE clover_min

SUBROUTINE clover_check_error(error)

  IMPLICIT NONE

  INTEGER :: error

  INTEGER :: maximum

  INTEGER :: err

  maximum=error

  CALL MPI_ALLREDUCE(error,maximum,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,err)

  error=maximum

END SUBROUTINE clover_check_error


END MODULE clover_module
