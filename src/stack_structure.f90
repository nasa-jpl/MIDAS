MODULE stack_structure
    IMPLICIT NONE
    Integer, PARAMETER :: max_size = 10000
    REAL(8),DIMENSION(max_size) :: stack
    INTEGER :: stack_top =0 ;
    LOGICAL :: stack_full, stack_empty
    
    CONTAINS
        SUBROUTINE push(pushed_no)
            REAL(8), INTENT(IN) :: pushed_no
            LOGICAL :: stack_full, stack_empty
            If (stack_full) Then
                print *, 'stack full. cannot push no'
                return
            EndIf
            stack_top = stack_top + 1
            stack(stack_top) = pushed_no
            stack_empty = .false.
            If (stack_top == max_size) Then
                stack_full= .true.
            EndIf           
        END SUBROUTINE push
        
        SUBROUTINE pop(popped_no)
            REAL(8), INTENT(OUT) :: popped_no
            If (stack_empty) Then
                print *, 'stack empty. cannot pop no'
                return
            EndIf
            popped_no = stack(stack_top)
            stack_top = stack_top - 1
            stack_full= .false.
            If (stack_top < 0) Then
                stack_empty= .true.
            EndIf           
        END SUBROUTINE pop
        
END MODULE stack_structure