package hageldave.optisled.generic.numerics;

import java.util.Arrays;
import java.util.Random;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.IntStream;

/**
 * Interface for operations on a generic matrix type {@code M}.
 * This allows for the implementation of numeric algorithms
 * (e.g. gradient descent) that are agnostic to the choice 
 * of a linear algebra library, and makes it possible to 
 * integrate this API into environments where a matrix
 * type is already present.
 * <p></p>
 * This collection of methods is by no means complete, but it contains
 * some of the more frequently used routines in linear algebra besides the basic operations,
 * such as EVD, SVD, Cholesky, pairwise distances, and some more.
 * When needing to express a loss function for a specific problem and need other operations,
 * you can extend this interface to add the needed methods
 * or implement the loss with a baked in linear algebra library
 * (in which case you loose the ability to use a drop-in replacement of MatCalc).
 * 
 * @author hageldave
 * @param <M> the matrix type
 */
 public interface MatCalc<M> {
	
	/**
	 * @param v values
	 * @return column vector of given values
	 */
	 M vecOf(double ... v);
	
	/**
	 * @param nRows number of rows
	 * @param values (in row major order)
	 * @return matrix with {@code nRows} and {@code nCols=values.length/nRows} filled with specified values
	 */
	 M matOf(int nRows, double... values);
	
	/**
	 * @param values  of the matrix
	 * @return matrix from the specified 2d array {@code (rowValues[] = values[i])}
	 */
	 M matOf(double[][] values);
	
	/**
	 * @param n number of dimensions
	 * @return zero column vector of dimension {@code n}
	 */
	 M zeros(int n);
	
	/** 
	 * @param rows number of rows
	 * @param columns number of columns
	 * @return matrix of zeros of specified size
	 */
	 M zeros(int rows, int columns);
	
	/**
	 * @param n number of rows (=number of columns)
	 * @return identity matrix of size n x n
	 */
	 default M eye(int n) {
		return eye(n, 1.0);
	}
	
	/**
	 * @param n number of rows (=number of columns)
	 * @param s scaling factor
	 * @return s*identity matrix of size n x n 
	 */
	 M eye(int n, double s);
	
	/**
	 * @param rows number of rows
	 * @param cols number of columns
	 * @param rnd random number generator
	 * @return matrix of random values in range 0.0 .. 1.0
	 */
	 default M rand(int rows, int cols, Random rnd) {
		M m = zeros(rows, cols);
		for(int r=0; r<rows; r++)
			for(int c=0; c<cols; c++)
				set_inp(m, r, c, rnd.nextDouble());
		return m;
	}
	
	/** 
	 * @param rows number of rows
	 * @param cols number of columns
	 * @return matrix of random values in range 0.0 .. 1.0
	 */
	 default M rand(int rows, int cols) {
		return rand(rows, cols, new Random());
	}
	
	/**
	 * @param n
	 * @return column vector of random values in range 0.0 .. 1.0
	 */
	 default M rand(int n) {
		return rand(n, 1);
	}
	
	/**
	 * @param rows number of rows
	 * @param cols number of columns
	 * @param rnd random number generator
	 * @return matrix with normally distributed random values (standard normal, mean=0, var=1)
	 */
	 default M randN(int rows, int cols, Random rnd) {
		M m = zeros(rows, cols);
		for(int r=0; r<rows; r++)
			for(int c=0; c<cols; c++)
				set_inp(m, r, c, rnd.nextGaussian());
		return m;
	}
	
	/**
	 * @param rows number of rows
	 * @param cols number of columns
	 * @return random matrix with normally distributed values (standard normal, mean=0, var=1)
	 */
	 default M randN(int rows, int cols) {
		return randN(rows, cols, new Random());
	}
	
	/**
	 * @param m matrix
	 * @return number of rows
	 */
	 int numRows(M m);
	
	/**
	 * @param m matrix 
	 * @return number of columns */
	 int numCols(M m);
	
	/**
	 * @param m matrix
	 * @return number of elements
	 */
	 default int numElem(M m) {
		return numCols(m)*numRows(m);
	}
	
	/**
	 * @param a vector 
	 * @param b vector
	 * @return inner or dot product for given vectors ‹a,b›
	 */
	 double inner(M a, M b);
	
	/**
	 * just an alias for {@link #inner(M, M)}
	 * @param a vector 
	 * @param b vector
	 * @return inner or dot product for given vectors ‹a,b›
	 */
	 default double dot(M a, M b) {return inner(a,b);}
	
	/**
	 * @param m matrix/vector
	 * @param s scaling factor
	 * @return a scaled version of given matrix/vector
	 */
	 M scale(M m, double s);
	
	/**
	 * @param m matrix/vector to be scaled in-place
	 * @param s scaling
	 * @return the same matrix/vector argument which was scaled in-place 
	 */
	 M scale_inp(M m, double s);
	
	/**
	 * @param a matrix/vector
	 * @param b matrix/vector
	 * @return result of matrix multiplication {@code a * b}
	 */
	 M matmul(M a, M b);
	
	/**
	 * @param a matrix/vector
	 * @param b matrix/vector
	 * @return element-wise multiplication of a and b (need to be of same size)
	 */
	 M elmmul(M a, M b);
	
	/**
	 * @param m matrix
	 * @param colV column vector
	 * @return copy of m where each row i was scaled by entry i of column vector
	 */
	 M mulRowsByColVec(M m, M colV);
	
	/**
	 * @param m matrix
	 * @param rowV row vector
	 * @return copy of m where each column j was scaled by entry j of row vector
	 */
	 M mulColsByRowVec(M m, M rowV);
	
	/**
	 * @param m matrix
	 * @return column vector of row sums of m
	 */
	 M rowSums(M m);
	
	/**
	 * @param m matrix
	 * @return column vector of row means of m
	 */
	 default M rowMeans(M m) {
		return scale_inp(rowSums(m), 1.0/numCols(m));
	}
	
	/**
	 * @param m matrix
	 * @return row vector of column means of m
	 */
	 default M colMeans(M m) {
		return scale_inp(colSums(m), 1.0/numRows(m));
	}
	
	/**
	 * @param m matrix
	 * @return column vector of row minimums
	 */
	 M rowMins(M m);
	
	/**
	 * @param m matrix
	 * @return column vector of row maximums
	 */
	 M rowMaxs(M m);
	
	/**
	 * @param m matrix
	 * @return row vector of column minimums
	 */
	 M colMins(M m);
	
	/**
	 * @param m matrix
	 * @return row vector of column maximums
	 */
	 M colMaxs(M m);
	
	/**
	 * @param m matrix
	 * @return row vector of column sums of m
	 */
	 M colSums(M m);
	
	/**
	 * @param m matrix
	 * @param rowV row vector
	 * @return result of subtracting specified row vector from each row of m
	 */
	 M subRowVec(M m, M rowV);
	
	/**
	 * @param m matrix
	 * @param colV column vector
	 * @return result of subtracting specified column vector from each column of m
	 */
	 M subColVec(M m, M colV);
	
	/**
	 * @param m matrix
	 * @param rowV row vector
	 * @return result of adding specified row vector to each row of m
	 */
	 M addRowVec(M m, M rowV);
	
	/**
	 * @param m matrix
	 * @param colV column vector
	 * @return result of adding specified column vector to each column of m
	 */
	 M addColVec(M m, M colV);
	
	/**
	 * @param m matrix
	 * @return transpose of m
	 */
	 M trp(M m);
	
	/**
	 * @param a matrix/vector
	 * @param b scalar
	 * @return result of adding b to each entry of a
	 */
	 M add(M a, double b);
	
	/**
	 * @param a matrix/vector
	 * @param b matrix/vector
	 * @return a+b (need to be same size)
	 */
	 M add(M a, M b);
	
	/**
	 * @param a matrix/vector
	 * @param b matrix/vector
	 * @return a-b (in-place subtraction from a, need to be same size)
	 */
	 M sub(M a, M b);
	
	/**
	 * @param a matrix/vector
	 * @param b matrix/vector
	 * @return a=a+b (in-place addition to a, need to be same size)
	 */
	 M add_inp(M a, M b);
	
	/**
	 * @param a matrix/vector
	 * @param b matrix/vector
	 * @return a=a-b (in-place subtraction from a, need to be same size)
	 */
	 M sub_inp(M a, M b);
	
	/**
	 * @param a matrix/vector
	 * @param b scalar
	 * @return result of subtracting b from each entry of a
	 */
	 default M sub(M a, double b) {
		return add(a, -b);
	}

	/**
	 * @param m matrix/vector
	 * @param idx index
	 * @param v value
	 * @return the same matrix/vector argument where the entry at given index was set (row major indexing)
	 */
	 M set_inp(M m, int idx, double v);

	/**
	 * @param m matrix/vector
	 * @param row row index
	 * @param col column index
	 * @param v value
	 * @return the same matrix/vector argument where the entry at specified position was set
	 */
	 default M set_inp(M m, int row, int col, double v) {
		return set_inp(m, row*numCols(m)+col, v);
	}

	/**
	 * @param m matrix/vector
	 * @param idx index
	 * @return entry of matrix/vector at given index (row major indexing)
	 */
	 double get(M m, int idx);

	/**
	 * @param m matrix/vector
	 * @param row row index
	 * @param col column index
	 * @return entry of matrix/vector at specified position
	 */
	 default double get(M m, int row, int col) {
		return get(m, row*numCols(m)+col);
	}

	/**
	 * @param m matrix/vector
	 * @return a copy of the argument
	 */
	 M copy(M m);
	
	/**
	 * @param m matrix/vector
	 * @return values of the matrix/vector in row major order */
	 double[] toArray(M m);
	
	/**
	 * @param m matrix/vector
	 * @return values of the matrix/vector {@code ( values[row][col] )} */
	 double[][] toArray2D(M m);
	
	/**
	 * @param m vector
	 * @param thresh norm threshold for when a vector is considered to be 0 (and thus not normalized)
	 * @return same as argument which was normalized in-place to unit length (if thresh smaller than norm) */
	 default M normalize_inp(M m, double thresh) {
		double norm = norm(m);
		return norm < thresh ? m : scale_inp(m, 1/norm);
	}
	
	/**
	 * @param m vector
	 * @param thresh norm threshold for when a vector is considered to be 0 (and thus not normalized)
	 * @return normalized to unit length version of vector (if thresh smaller than norm) */
	 default M normalize(M m, double thresh) { return normalize_inp(copy(m), thresh); }
	
	/**
	 * @param m vector
	 * @return same as argument which was normalized in place to unit length (if 1e-7 smaller than norm) */
	 default M normalize_inp(M m) { return normalize_inp(m, 1e-7); }
	
	/**
	 * @param m vector
	 * @return normalized to unit length version of vector (if 1e-7 smaller than norm) */
	 default M normalize(M m) { return normalize_inp(copy(m)); }
	
	/**
	 * @param v vector
	 * @return squared vector norm ||v||^2 = ‹v,v› */
	 default double norm2(M v) { return inner(v,v); }
	
	/**
	 * @param v vector
	 * @return vector norm ||v|| = sqrt(‹v,v›) */
	 default double norm(M v) { return Math.sqrt(norm2(v)); }

	/**
	 * @param m matrix
	 * @return squared Frobenius norm of the matrix (sum of squared elements)
	 */
	 default double frob2(M m) {
		return sum(elmmul(m, m));
	}

	/**
	 * @param m matrix
	 * @return Frobenius norm of the matrix (sqrt of sum of squared elements)
	 */
	 default double frob(M m) {
		return Math.sqrt(frob2(m));
	}
	
	/**
	 * @param m matrix
	 * @param full when false, allows omitting superfluous components, otherwise U and V will be computed as complete bases
	 * @return singular value decomposition (full or sparse/economic) with matrices U,S,V where m=U*S*trp(V).
	 * singular values are ordered from largest to smallest along the diagonal (largest at 0,0)
	 */
	 M[] svd(M m, boolean full);
	
	/**
	 * @param m matrix
	 * @returns Cholesky decomposition of m, the upper triangular part U, so that m = trp(U)*U */
	 M cholesky(M m);

	/**
	 * @param m matrix (symmetric)
	 * @return eigendecomposition for symmetric matrices, with matrices Q,S where m = Q*S*Q^-1 (eigenvectors in columns of Q)
	 */
	 M[] symEvd(M m);
	
	/**
	 * @param m matrix
	 * @return determinant of m */
	 double det(M m);
	
	/**
	 * @param m matrix
	 * @return element-wise exp (in-place) */
	 default M exp_inp(M m) {
		return elemwise_inp(m, Math::exp);
	}
	
	/**
	 * @param m matrix
	 * @return element wise sqrt (in-place) */
	 default M sqrt_inp(M m) {
		return elemwise_inp(m, Math::sqrt);
	}
	
	/**
	 * @param m matrix
	 * @param f funtion to be applied to each element
	 * @return element wise f(m_ij) (in-place) */
	 default M elemwise_inp(M m, DoubleUnaryOperator f) {
		for(int i=0; i<numElem(m); i++)
			set_inp(m, i, f.applyAsDouble(get(m, i)));
		return m;
	}
	
	/**
	 * @param a matrix/vector
	 * @param b matrix/vector
	 * @return side by side concatenation (equal number of rows) */
	 M concatHorz(M a, M b);
	
	/**
	 * @param a matrix/vector
	 * @param b matrix/vector
	 * @return stacked concatenation (equal number of columns) */
	 M concatVert(M a, M b);

	/**
	 * @param m matrix/vector
	 * @param ra row start
	 * @param rb row end (exclusive)
	 * @param ca column start
	 * @param cb column end (exclusive)
	 * @return sub-matrix of rows ra to rb (exclusive) and columns ca to cb (exclusive)
	 */
	 M getRange(M m, int ra, int rb, int ca, int cb);

	/**
	 * @param m matrix/vector
	 * @param r row index
	 * @return row of m
	 */
	 default M getRow(M m, int r) {
		return getRange(m, r, r+1, 0, numCols(m));
	}

	/**
	 * @param m matrix/vector
	 * @param c column index
	 * @return column of m
	 */
	 default M getCol(M m, int c) {
		return getRange(m, 0, numRows(m), c, c+1);
	}
	
	/**
	 * @param m matrix
	 * @return the diagonal of the matrix as column vector */
	 M diagV(M m);
	
	/**
	 * @param v vector
	 * @return the diagonal matrix from the vector (vector elements put on diagonal) */
	 M diagM(M v);

	/**
	 * copies values of src to target (copies values in row-major order, matrices don't need to be of same size)
	 * @param src source matrix
	 * @param target target matrix
	 */
	 default void copyValues(M src, M target) {
		copyValues(src, 0, target, 0, Math.min(numElem(src), numElem(target)));
	}

	/**
	 * copies specified range of values from src to target
	 * (copies values in row-major order, matrices don't need to be of same size)
	 * @param src source matrix
	 * @param startSrc start index in source
	 * @param target target matrix
	 * @param startTarget start index in target
	 * @param len number of values to copy (row-major order)
	 */
	 default void copyValues(M src, int startSrc, M target, int startTarget, int len) {
		for(int i=0; i<len; i++)
			set_inp(target, i+startTarget, get(src, i+startSrc));
	}

	/**
	 * @param n array length
	 * @return array of type M (since Java does not allow creation of generically typed arrays)
	 */
	 M[] matArray(int n);

	/**
	 * @param m outer array length
	 * @param n inner array length
	 * @return 2D array of type M (since Java does not allow creation of generically typed arrays)
	 */
	 M[][] matArray(int m, int n);

	/**
	 * @param a matrix/vector
	 * @param b matrix/vector
	 * @return matrix multiplication y=a*b
	 */
	 default M mult_ab(M a, M b) {
		return matmul(a, b);
	}

	/**
	 * @param a matrix/vector
	 * @param b matrix/vector
	 * @return matrix multiplication y=trp(a)*b
	 */
	 M mult_aTb(M a, M b);

	/**
	 * @param a matrix/vector
	 * @param b matrix/vector
	 * @return matrix multiplication y=a*trp(b)
	 */
	 M mult_abT(M a, M b);

	/**
	 * @param a matrix/vector
	 * @param b matrix/vector
	 * @param c matrix/vector
	 * @return matrix multiplication y=trp(a)*b*c
	 */
	 default M mult_aTbc(M a, M b, M c) {
		return matmul(mult_aTb(a, b), c);
	}

	/**
	 * @param a matrix/vector
	 * @param b matrix/vector
	 * @param c matrix/vector
	 * @return matrix multiplication y=a*b*trp(c)
	 */
	 default M mult_abcT(M a, M b, M c) {
		return matmul(a, mult_abT(b, c));
	}

	/**
	 * @param a matrix/vector
	 * @param b matrix/vector
	 * @param c matrix/vector
	 * @return matrix multiplication y=a*trp(b)*c
	 */
	 default M mult_abTc(M a, M b, M c) {
		return matmul(a, mult_aTb(b, c));
	}

	/**
	 * @param m matrix/vector
	 * @return sum of elements
	 */
	 double sum(M m);

	/**
	 * @param v1 vector
	 * @param v2 vector
	 * @return squared distance d=||v1-v2||^2
	 */
	 default double dist2(M v1, M v2) {
		double sum=0;
		for(int i=0; i<numElem(v1); i++) {
			double diff=get(v1, i)-get(v2, i);
			sum += diff*diff;
		}
		return sum;
	}

	/**
	 * @param v1 vector
	 * @param v2 vector
	 * @return distance d=||v1-v2||
	 */
	default double dist(M v1, M v2) {
		return Math.sqrt(dist2(v1,v2));
	}

	/**
	 * @param a matrix of row vectors
	 * @param b matrix of row vectors
	 * @return pairwise distances between rows of a and b.
	 * Entry at (r,c) is distance between a_c and b_r.
	 */
	 default M pairwiseDistances2(M a, M b) {
		M dists = zeros(numRows(b), numRows(a));
		for(int i=0; i<numRows(b); i++) {
			M bi = getRow(b, i);
			for(int j=0; j<numRows(a); j++) {
				M aj = getRow(a, j);
				double distSquared = dist2(bi, aj);
				set_inp(dists, i, j, distSquared);
			}
		}
		return dists;
	}

	/**
	 * @param evd eigendecomposition [Q,S]
	 * @return sorted eigenvectors and eigenvalues in eigenvalue ascending order (in-place sorting)
	 */
	default M[] sortEVD_inp(M[] evd) {
		double[] negvals = toArray(scale(diagV(evd[1]), -1));
		int[] order = argsort(negvals);
		M vectors = copy(evd[0]);
		for(int i=0; i<order.length; i++) {
			int j = order[i];
			M ev = getCol(vectors, j);
			for(int r=0; r<numElem(ev); r++) {
				set_inp(evd[0], r, i, get(ev, r));
			}
			set_inp(evd[1], i, i, -negvals[j]);
		}
		return evd;
	}

	/**
	 * argument sorting
	 * @param toSort values to sort (will stay untouched)
	 * @return order array of indices. smallest=toSort[order[0]], largest=toSort[order[toSort.length-1]]
	 */
	 public static int[] argsort(final double[] toSort) {
		Integer[] indices = IntStream.range(0, toSort.length).mapToObj(Integer::valueOf).toArray(Integer[]::new);
		Arrays.sort(indices, (i,j)->Double.compare(toSort[i], toSort[j]));
		return Arrays.stream(indices).mapToInt(Integer::intValue).toArray();
	}

	/**
	 * @param m matrix
	 * @return pseudo inverse of m
	 */
	 default M pinv(M m) {
		M[] svd = svd(m, false);
		M Sinv = diagM(elemwise_inp(diagV(svd[1]), v->1.0/v));
		return mult_abcT(svd[2],Sinv,svd[0]);
	}

	/**
	 * converts a matrix of this {@link MatCalc} to a matrix for use with another MatCalc, e.g.
	 * when using multiple linear algebra libraries.
	 * @param m matrix
	 * @param otherMC other MatrixCalculator
	 * @param <Other> other matrix type
	 * @return a copy of matrix m but as type of otherMC.
	 */
	 default <Other> Other convert(M m, MatCalc<Other> otherMC){
		int rows = numRows(m), cols = numCols(m);
		Other o = otherMC.zeros(rows,cols);
		for(int i=0; i<rows*cols; i++){
			otherMC.set_inp(o,i,get(m,i));
		}
		return o;
	}
}










