with(LinearAlgebra):

# This is just an abbreviation:
RREF := ReducedRowEchelonForm:
RO := RowOperation:
CO := ColumnOperation:

# The symbol D is used by default for differentiation, so we need to tell Maple
# that we will be using it as the ame of a matrix:
unprotect(D);

# This function checks whether a matrix is in Reduced Row Echelon Form
is_rref := proc(A) Equal(RREF(A),A); end:

# This function converts an augmented matrix to a system of equations
augmat_to_eqs := proc(A,vars)
 local L,R;
 L := A . Vector([op(vars),0]);
 R := A . Vector([ 0 $ (ColumnDimension(A) - 1),1]);
 [seq(L[i] = R[i],i=1..RowDimension(A))];
end:

# This function converts a system of equations to an augmented matrix
eq_to_row := proc(e,vars)
 local p,x;
 p := lhs(e) - rhs(e);
 return [seq(coeff(p,x,1),x in vars),-subs({seq(x = 0,x in vars)},p)];
end:

eqs_to_augmat := proc(E,vars)
 return(Matrix(map(eq_to_row,[op(E)],vars)));
end:

# This function takes two vectors u and v, and produces a list of equations saying that
# the components of u and v are equal.
component_eqs := proc(u,v)
 [seq(u[i]=v[i],i=1..Dimension(u))];
end:

# This function carries out a sequence of row operations, showing all the intermediate 
# matrices that are produced along the way.
show_rowops := proc(A) 
 local B,x,y,hide;
 B := Copy(A);
 print(B);
 for x in args[2..-1] do
  if type(x,function) and op(0,x) = H then
   hide := true;
  else
   hide := false;
  fi;
  
  RowOperation(B,op(x),inplace=true);
  B := map(expand,B);
  if not(hide) then print(B); fi;
 od;
 return(eval(B));
end:

# In MAS201 we define the characteristic polynomial to be the determinant of A-tI.
# The built-in Maple function gives the determinant of tI-A instead.  Here we define
# a function that follows the MAS201 convention.
char_poly := (A,t) -> 
  sort((-1)^RowDimension(A) * CharacteristicPolynomial(A,t)):

# For convenience, we give names to the identity matrices of various sizes.
I2 := IdentityMatrix(2):
I3 := IdentityMatrix(3):
I4 := IdentityMatrix(4):
I5 := IdentityMatrix(5):

# This function rotates a matrix or vector through 180 degrees, which is one ingredient
# in our method for computing the canonical basis for the annihilator of a list of vectors.
rotate := proc(A)
 local n,m,i,j,B;
 if type(A,Matrix) then
  n := RowDimension(A);
  m := ColumnDimension(A);
  B := Matrix(n,m);
  for i from 1 to n do
   for j from 1 to m do
    B[i,j] := A[n+1-i,m+1-j];
   od;
  od;
  return(B);
 elif type(A,Vector) then
  n := Dimension(A);
  B := Vector(n);
  for i from 1 to n do
   B[i] := A[n+1-i];
  od;
  return(B);
 fi;
end:

# This function is used to help sort a list of vectors into canonical order.
compare_vectors := proc(a,b)
 local i;
 for i from 1 to Dimension(a) do
  if (abs(a[i]) <> abs(b[i])) then return(evalb(abs(a[i])>abs(b[i]))); fi;
  if a[i] <> b[i] then return(evalb(a[i]>b[i])); fi;
 od;
 return(false);
end:

# This function computes the canonical basis for the span of a list of vectors. 
span_basis := proc()
 ColumnSpace(Matrix([args]));
end:

# This function computes the canonical basis for the annihilator of a list of vectors
ann_basis := proc()
 local U,L;
 U := Matrix([args]);
 L := [op(map(rotate,NullSpace(rotate(Transpose(U)))))];
 L := sort(L,compare_vectors);
 return(L);
end:

# The various elmat functions produce elementary matrices.  They all assume that the
# global variable n has been set to the desired size of the matrices.
elmat_D := proc(p,lambda)
 global n;
 local A,i;
 A := Matrix(n,n);
 for i from 1 to n do A[i,i] := 1; od;
 A[p,p] := lambda;
 return(eval(A));
end:

elmat_E := proc(p,q,mu)
 global n;
 local A,i;
 if p=q then
  return(FAIL);
 fi;
 A := Matrix(n,n);
 for i from 1 to n do A[i,i] := 1; od;
 A[p,q] := mu;
 return(eval(A));
end:

elmat_F := proc(p,q)
 global n;
 local A,i;
 if p=q then
  return(FAIL);
 fi;
 A := Matrix(n,n);
 for i from 1 to n do A[i,i] := 1; od;
 A[p,p] := 0;
 A[q,q] := 0;
 A[p,q] := 1;
 A[q,p] := 1;
 return(eval(A));
end:

# The following function computes the matrix that we have called the adjugate.
# Some authors call it the 'classical adjoint' instead, and Maple's naming reflects this.
adjugate := Adjoint:

# The following function constructs diagonal matrices.
diag := () -> DiagonalMatrix([args]):

starting_slot := proc(u) 
 local n,i;
 if type(u,list) then
  n := nops(u);
 else
  n := Dimension(u);
 fi;
 for i from 1 to n do
  if u[i] <> 0 then
   return(i);
  fi;
 od;
 return(infinity);
end:

matrix_latex := proc(A)
 local n,m,i,j,s;
 n := RowDimension(A);
 m := ColumnDimension(A);
 s := "\\bbm\n";
 for i from 1 to n do
  for j from 1 to m do
   if j>1 then 
    s := cat(s,"&",sprintf("%Q",A[i,j]));
   else 
    s := cat(s,sprintf("%Q",A[i,j]));
   fi;
  od;
  s := cat(s,"\\\\\n");
 od;
 s := cat(s,"\\ebm\n");
 return(s);
end: 

show_rowops_latex := proc(A) 
 local B,x,s,hide;
 s := "\\[\n";
 B := Copy(A);
 B := map(expand,B);
 s := cat(s,matrix_latex(B));
 for x in args[2..-1] do
  if type(x,function) and op(0,x) = H then
   hide := true;
  else
   hide := false;
  fi;

  RowOperation(B,op(x),inplace=true);
  B := map(expand,B);
  if not(hide) then s := cat(s,"\\to\n",matrix_latex(B)); fi;
 od;
 s := cat(s,"\n\\]\n");
 printf(s);
end:

nice_invertible_matrix := proc(n,d)
 local A,r,i,j,k;
 r := rand(-d..d);
 A := Matrix(n,n); 
 k := 100;

 while(Determinant(A)^2 <> 1 and k > 0) do 
  A := Matrix(n,n,[seq(r(),i=1..n^2)]);
  k := k - 1;
 od;

 if k = 0 then
  return(FAIL);
 else 
  return(A);
 fi;
end:
