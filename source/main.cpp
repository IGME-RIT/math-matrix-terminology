/*
Title: Matrix Mathematics
File Name: main.cpp
Copyright © 2016
Author: Andrew Litfin
Written under the supervision of David I. Schwartz, Ph.D., and
supported by a professional development seed grant from the B. Thomas
Golisano College of Computing & Information Sciences
(https://www.rit.edu/gccis) at the Rochester Institute of Technology.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.
This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// The primary objects of study in linear algebra are matrices.
// This tutorial series will explore the applications of matrices to computer games and simulation,
//  especially in the realm of physical transformations.
// The exposition follows that of Eric Lengyel in "Foundations of Game Engine Development" (Volume 1).
// We have included the Vector structs from the previous series, and introduced Matrix structs that act similarly.
// These structs are based upon and largely follow code samples given in FGED.
//  As before, Matrix2D is heavily annotated, with other structs being annotated in places of difference.

// In this tutorial, we examine some basic terminology regarding matrices.

#include "../header/Matrix4D.h"
#include "../header/tests.h"
#include "../header/helpers.h"

#include <iostream>
#include <ctime>

int main()
{
	// Required for some helper functions
	srand((unsigned)time(0));

	// Square matrices
	// ---------------
	// We call a matrix "square" when it has the same number of rows as columns.
	// Every matrix in these examples is square, but you should know that these are not the general case.

	// Diagonal elements
	// -----------------
	// The diagonal elements of a matrix are (shockingly) the elements on the main diagonal of the matrix,
	//  that is, the elements A(i, i).
	// Note that a matrix need not be square to refer to the "diagonal" elements.
	Matrix3D A(0, 1, 2, 3, 4, 5, 6, 7, 8);
	std::cout << "A =\n" << A << "0, 4, and 8 are the only diagonal elements of A\n";

	// Diagonal matrix
	// ---------------
	// A diagonal matrix is a matrix with 0 on every element of the matrix that is not a diagonal element.
	Matrix3D diagonal(1, 0, 0, 0, 2, 0, 0, 0, 3);
	std::cout << "A diagonal matrix:\n" << diagonal;

	// Transpose
	// ---------
	// The Transpose of a matrix is the matrix achieved by swapping the rows and columns.
	// Technically, given a matrix A, then Transpose(A)(i, j) = A(j, i), by definition.
	// The transpose is most often denoted with a superscript T, but that is difficult to achieve in a monospace font.
	std::cout << "Transpose(A) =\n" << Transpose(A);

	// Symmetric
	// ---------
	// A matrix is symmetric if A(i, j) = A(j, i), or equivalently, if A = Transpose(A).
	// Note that a matrix must be square to be symmetric, and that every diagonal matrix is symmetric automatically.
	Matrix3D symmetric(1, 2, 3, 2, 4, 5, 3, 5, 6);
	std::cout << "A symmetric matrix:\n" << symmetric;

	// Skew-symmetric
	// --------------
	// A matrix is skew-symmetric if its transpose is its negation, i.e. Transpose(A) = -A.
	// As an immediate consequence, all diagonal elements are 0 (being that 0 is the only number equal to its negative).
	// Skew-symmetric matrices are a fascinating subject in general linear algebra, but for the purposes of game design they
	//  are almost never used for purposes other than representing the cross product as a matrix.
	// That is,
	//  Cross(a, b) = CrossMat(a) * b,
	//  where CrossMat(a) is a skew-symmetric matrix.
	Vector3D a(randFloat(-1, 1), randFloat(-1, 1), randFloat(-1, 1));	
	Vector3D b(randFloat(-1, 1), randFloat(-1, 1), randFloat(-1, 1));
	std::cout << "a = " << a << ", b = " << b << "\n"
		"Cross(a, b) = " << Cross(a, b) << "\n"
		"CrossMat(a) =\n" << CrossMat(a) <<
		"CrossMat(a) * b = " << CrossMat(a) * b << "\n";
	if (Cross(a, b) == CrossMat(a) * b)
	{
		std::cout << "They are equivalent!\n";
	}

	// Identity Matrix
	// ---------------
	// The matrix in which every diagonal element is 1 and every off diagonal element is 0.
	// This matrix has the unique property that for all matrices A, I*A = A*I = A.
	// It is called the "identity matrix" because it is the multiplicative identity.
	std::cout << "The identity matrix:\n" << Matrix3D();

	// Orthogonal
	// ----------
	// A matrix is orthogonal if its transpose is its inverse, i.e. Q*Transpose(Q) = Transpose(Q)*Q = I
	// By definition, then, det(Q)^2 = 1, hence det(Q) = +1 or -1.
	// An example of an orthogonal matrix is rotation matrices.
	Matrix3D rotZ = MakeRotationZ(3.1415f / 4);
	std::cout << "rotX =\n" << rotZ
		<< "rotX*Transpose(rotZ) =\n" << rotZ * Transpose(rotZ)
		<< "Transpose(rotZ) =\n" << Transpose(rotZ) * rotZ
		<< "Det(rotZ) = " << Determinant(rotZ) << "\n";

	// Idempotent
	// ----------
	// A matrix is idempotent if its square is equal to itself, i.e. P*P = P.
	// Necessarily, only square matrices can be idempotent.
	// This is in fact the defining property of a projection matrix.
	Matrix3D P = MakeProjection(Vector3D(1, 0, 0)) + MakeProjection(Vector3D(0, 1, 0));
	std::cout << "P =\n" << P
		<< "P*P =\n" << P*P;
	if (P*P == P)
	{
		std::cout << "P is idempotent!\n";
	}

	// Involution
	// ----------
	// A matrix is an involution if its square is the identity, i.e. T*T = I.
	// A common involution is a reflection.
	Matrix3D T = Matrix3D() - 2 * Outer(Vector3D(1, 0, 0), Vector3D(1, 0, 0));
	std::cout << "T =\n" << T
		<< "T*T =\n" << T*T;
	if (T*T == Matrix3D())
	{
		std::cout << "T is an involution!\n";
	}

	// Image
	// -----
	// The image of any function f : A -> B is defined as f[A' \subseteq A] = { f(a) | a \in A' }.
	// Informally, the image of a function is the set of possible values taken by a function.
	// Since matrices are a specific type of function, this definition carries over.
	// More often though, the image of a matrix is notated as im(A).
	// The image of a matrix is itself a vector space, in that A(su + v) = sA(u) + A(v)
	// In particular, there is an obvious choice of basis for that vector space:
	//  the columns of A.
	// Then im(A) = span({ A[0], A[1], ... A[n-1] }) = { Sum_{i}(x_i * A[i]) | x_i \in R }

	// Determinant
	// -----------
	{
		/*
		The determinant of a square matrix over the reals is a real number that is the (signed) volume of the parallelotope
		spanned by the columns of the matrix.
		It is a measure of how much a matrix stretches or compresses space.
		If it is positive, then the parallelotope has the same orientation as the surrounding space.
		If it is negative, then the parallelotope has the reverse orientation.
		(By orientation, we mean "left-handed" versus "right-handed" space, informally.)

		Matrices with a determinant of 0 are called "singular" matrices, and are rather annoying.
		A matrix has an inverse if and only if it is nonsingular.

		There are many ways to calculate the determinant.
		The highest level (and least often used) is called Leibniz expansion.
		To calculate the determinant this way, consider the set of permutations on n elements, say {0, 1, 2, ..., n-1}.
		There are n! such permutations. Then sum over all such permutations, taking the product of the elements of the matrix,
		 where the row element is an index and the column element is the resulting element given by the permutation.
		The set of all such permutations is called the symmetric group on n elements and is denoted S_n.
		Then the Leibniz expansion is

		det(A) = Sum_{ s \in S_n }(sgn(s) * Prod_{i=0}^{n-1}(A(i, s(i))))

		where sgn denotes the sign of the permutation.
		sgn(s) = +1 if the permutation is even, and -1 if the permutation is odd.

		This is not particularly good for computation, so instead there is an equivalent method known as Minor expansion.
		A matrix minor of an n by n matrix is an (n-1) by (n-1) matrix excluding a given row and column.
		It is a recursive definition as follows:

		Choose a fixed k. Then
		det(A) = Sum_{j=0}^{n-1}(A(k, j) * (-1)^{k + j} * det(Minor(A, k, j)))

		For small matrices, you will most often see it expanded.
		E.G., for a 2D matrix,
		A =
		[ a b ]
		[ c d ]
		det(A) = a*d - b*c

		The determinant has many nice properties:
		  det(I) = 1, where I is the identity matrix
		  det(Transpose(A)) = det(A)
		  det(AB) = det(A)*det(B)
		  det(Inverse(A)) = 1/det(A)
		  det(t*A) = t^n * det(A), where t is any scalar and n is the dimension of the matrix
		*/
		std::cout << "Det(A) = " << Determinant(A) << "\n";
	}

	// Elementary Matrices
	// -------------------
	// Elementary matrices are called such because they are `atomic` in the sense that all other matrices can be realized as products of elementary matrices.
	// Note that the described effects are when left-multiplied, i.e. A' has the described effect when A' = E*A
	// There are 3 types of elementary matrices:
	// 1) Scaling matrices
	//    These are nearly identical to the identity matrix, but one diagonal element is neither 0 nor 1.
	std::cout << "Scale matrix:\n" << Matrix3D(2, 0, 0, 0, 1, 0, 0, 0, 1);
	// 2) Swap matrices
	//    As the name implies, these are identity matrices with two rows swapped.
	std::cout << "Swap matrix:\n" << Matrix3D(1, 0, 0, 0, 0, 1, 0, 1, 0);
	//    Swaps the second and third row of a matrix.
	// 3) Addition matrices
	//    Less obviously, these matrices add a scalar multiple of one row to another row.
	//    Specifically, when E(r, s) = t, it adds t*A.row(s) to A.row(r)
	std::cout << "Addition matrix:\n" << Matrix3D(1, 1, 0, 0, 1, 0, 0, 0, 1);

	// Inverse
	// -------
	{
		// The inverse of a matrix A is the unique matrix A^{-1} or Inverse(A) such that A*A^{-1} = A^{-1}*A = I
		// There are many ways to compute the inverse of a matrix.
		// For our purposes, we use algorithms tailored to each dimension of matrix, rather than a general case.

		// For example, given a 2D matrix
		// [ a b ]
		// [ c d ]
		// The inverse is
		//  __1__ [ d -b ]
		//  ad-bc [ -c a ]

		// However, there does exist a `general' case, which is included with these examples.
		// We need some additional definitions, however:

		Matrix3D m(0, 1, 2, -3, 4, 5, 6, 7, 8);
		std::cout << "m =\n" << m;

		// Minor
		// -----
		// The (i, j) minor of an n by n matrix is an n-1 by n-1 matrix excluding row i and column j.
		std::cout << "The (0, 0) minor of m is\n" << Minor(m, 0, 0);

		// Cofactor
		// --------
		// The (i, j) cofactor of a matrix is the determinant of the (i, j) minor times -1 if i+j is odd
		std::cout << "The (0, 0) cofactor of m is " << Cofactor(m, 0, 0) << ".\n";

		// Cofactor Matrix
		// ---------------
		// The cofactor matrix of a matrix is the matrix C where C(i, j) = Cofactor(m, i, j)
		std::cout << "The cofactor matrix of m is\n" << CofactorMatrix(m);

		// Adjugate
		// --------
		// The adjugate matrix of a matrix is simply the transpose of the cofactor matrix.
		std::cout << "The adjugate of m is\n" << Adjugate(m);

		// Then the Adjugate inverse formula is
		// InverseAdj(m) = Adjugate(m)/Determinant(m)
		std::cout << "InverseAdj(m) =\n" << InverseAdj(m);

		// Compare to the inverse computed `explicitly':
		std::cout << "Inverse(m) =\n" << Inverse(m);

		if (Inverse(m) == InverseAdj(m))
		{
			std::cout << "The two are equivalent!\n";
		}

		// Note, however, that this method is computationally more intensive.
		TestMat3DInverseSpeed(10000);

		// This method becomes unweildy very fast, so there exist inverse algorithms, such as Gauss-Jordan, to find the inverse.
	}

	std::cout << "\nPress Enter to exit . . . ";
	std::cin.get();
	return 0;
}
