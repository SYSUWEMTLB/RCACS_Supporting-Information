      Subroutine SINTER(X,A,Y,B,M,N)
      SAVE
C     VERSION(05/30/91)
C
C-----------------------------------------------------------------------
C        A SPECIAL CASE OF INTERP ....NO EXTRAPOLATION BELOW AND ABOVE DATA
C        THIS ROUTINE LINEARLY INTERPOLATES AN ARRAY B
C        X(M) MUST BE DESCENDING
C        A(X) GIVEN FUNCTION
C        B(Y) FOUND BY LINEAR INTERPOLATION AND EXTRAPOLATION
C        Y(N) THE DESIRED DEPTHS
C        M    THE NUMBER OF POINTS IN X AND A
C        N    THE NUMBER OF POINTS IN Y AND B
C-----------------------------------------------------------------------
C
      Dimension X(M), A(M), Y(N), B(N)
C
C-------- EXTRAPOLATION CASES ------------------------------------------
      Do 10 I = 1, N
        If (Y(I).GT.X(1)) B(I) = A(1) 
        If (Y(I).LT.X(M)) B(I) = A(M)
   10 Continue
C
C-------- INTERPOLATION CASES ------------------------------------------
      NM = M - 1
      Do 30 I = 1, N
        Do 20 J = 1, NM
          If (Y(I).LE.X(J).AND.Y(I).GE.X(J+1)) B(I) = A(J) - (A(J)-
     *        A(J+1)) * (X(J)-Y(I)) / (X(J)-X(J+1))
   20   Continue
   30 Continue
C
      Return
      End
