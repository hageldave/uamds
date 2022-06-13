package hageldave.optisled.ejml;

import java.util.Arrays;
import java.util.Locale;
import java.util.function.BiFunction;
import java.util.function.DoubleFunction;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import org.ejml.data.DMatrixRMaj;
import org.ejml.simple.SimpleMatrix;

public class MatUtil {

	/**
	 * Creates a new SimpleMatrix of specified dimensions, that is backed by
	 * the specified array (in row-major order). If the data array is null
	 * then a new array is allocated to back this matrix (all entries are zero then).
	 * 
	 * @param nrows number of rows
	 * @param ncols number of columns
	 * @param data optional (can be null) row-major order data array for the matrix
	 * @return nrows x ncols matrix
	 * 
	 * @throws IllegalArgumentException if nrows &lt; 1 or ncols &lt; 1 or when the 
	 * specified data array is not null but does not have nrows*ncols length.
	 */
	public static SimpleMatrix matrix(int nrows, int ncols, double... data){
		if(nrows < 1 )
			throw new IllegalArgumentException("number of rows has to be positive. nrows: " + nrows);
		if(ncols < 1 )
			throw new IllegalArgumentException("number of columns has to be positive. ncols: " + ncols);
		if(data == null)
			data = new double[nrows*ncols];
		else if(nrows*ncols != data.length)
			throw new IllegalArgumentException(String.format("the specified dimensions do not match the length of the data. [%dx%d]=%d elements, data array length = %d", nrows,ncols,ncols*nrows,data.length));
		
		return SimpleMatrix.wrap(DMatrixRMaj.wrap(nrows, ncols, data));
	}
	
	/**
	 * Creates a new SimpleMatrix of the specified dimensions where all entries are zero.
	 * 
	 * @param nrows number of rows
	 * @param ncols number of columns
	 * @return nrows x ncols matrix
	 * @throws IllegalArgumentException if nrows &lt; 1 or ncols &lt; 1
	 */
	public static SimpleMatrix matrix(int nrows, int ncols){
		return matrix(nrows, ncols, null);
	}
	
	/**
	 * Creates a new SimpleMatrix with 1 column and the specified number of rows
	 * (column vector). When the data array argument is not null it will be used
	 * to back the matrix, otherwise a new array will be allocated (all entries will
	 * be zero then).
	 * @param n number of elements (size of vector, number of rows)
	 * @param data optional data array for the vectors entries.
	 * @return n dimensional column vector
	 * @throws IllegalArgumentException if n &lt; 1 or when the specified data array
	 * is not null but has not length n.
	 */
	public static SimpleMatrix vector(int n, double[] data){
		return matrix(n, 1, data);
	}
	
	/**
	 * Creates a new column vector of the specified size, where all entries are zero.
	 * @param n number of elements (size of vector, number of rows)
	 * @return n dimensional column vector
	 * @throws IllegalArgumentException if n &lt; 1
	 */
	public static SimpleMatrix vector(int n){
		return vector(n, null);
	}
	
	/**
	 * Creates a new column vector of the specified size, where all entries have
	 * the specified value.
	 * @param n number of elements (size of vector, number of rows)
	 * @param value for entries of the vector
	 * @return n dimensional column vector with all entries of specified value
	 * @throws IllegalArgumentException if n &lt; 1
	 */
	public static SimpleMatrix vector(int n, double value){
		SimpleMatrix vector = vector(n);
		vector.fill(value);
		return vector;
	}
	
	/**
	 * Creates a new column vector of the specified size, where all entries are 1.
	 * @param n number of elements (size of vector, number of rows)
	 * @return n dimensional column vector with all entries 1
	 * @throws IllegalArgumentException if n &lt; 1
	 */
	public static SimpleMatrix ones(int n){
		return vector(n,1.0);
	}
	
	/**
	 * Applies the specified function f to each element of a copy of the matrix M.
	 * <br> M'<sub>ij</sub> = f(M<sub>ij</sub>), <i>(i=column, j=row)</i>
	 * 
	 * @param m the matrix 
	 * @param fn the function to apply
	 * @return a copy of m with fn applied to each element
	 */
	public static SimpleMatrix applyFn(SimpleMatrix m, DoubleFunction<Double> fn){
		m = m.copy();
		for(int i = 0; i < m.getNumElements(); i++){
			m.set(i, fn.apply(m.get(i)));
		}
		return m;
	}
	
	/**
	 * Applies the specified function f<sub>i</sub> to each column i of a copy of the matrix M.
	 * <br> M'<sub>ij</sub> = f<sub>i</sub>(M<sub>ij</sub>), <i>(i=column, j=row)</i> <br>
	 * 
	 * @param m the matrix
	 * @param fn the function to apply, the function inputs are fn(value, column).
	 * @return a copy of m with fn applied to each element of each column
	 */
	public static SimpleMatrix applyFnColumnWise(SimpleMatrix m, BiFunction<Double, Integer, Double> fn){
		m = m.copy();
		for(int col = 0; col < m.numCols(); col++){
			for(int row = 0; row < m.numRows(); row++){
				m.set(row, col, fn.apply(m.get(row, col), col));
			}
		}
		return m;
	}
	
	/**
	 * Adds to each element of each column of the matrix M the corresponding value of the vector s.
	 * <br> M'<sub>ij</sub> = M<sub>ij</sub> + s<sub>i</sub>, <i>(i=column, j=row)</i>
	 * @param m the matrix
	 * @param summand the row vector s (can also be a column vector)
	 * @return a copy of m with the ith element of the summand added to each element of the ith column of m
	 * @throws ArrayIndexOutOfBoundsException if the summand vector has less elements than the matrix has columns
	 */
	public static SimpleMatrix addColumnWise(SimpleMatrix m, SimpleMatrix summand){
		return applyFnColumnWise(m, (v,col)->v+summand.get(col));
	}
	
	/**
	 * Returns a String representation of the specified matrix's dimensions.
	 * e.g. {@code [3 rows x 2 cols]}
	 * @param m the matrix
	 * @return matrix dimensions as String
	 */
	public static String dimensionToString(SimpleMatrix m){
		return String.format("[%d rows x %d cols]", m.numRows(), m.numCols());
	}
	
	/** 
	 * @param cause the causing SingularMatrixException
	 * @param matrix the matrix which is singular or is subject to the cause of the exception
	 * @param message the message of the returned exception.
	 * @return a SingularMatrixException with attached matrix for debugging 
	 */
	public static SingularMatrixException singularMatrixException(org.ejml.data.SingularMatrixException cause, SimpleMatrix matrix, String message){
		return new SingularMatrixException(cause, matrix, message);
	}
	
	/** extension of SingularMatrixException that can store a matrix for debugging purposes */
	public static class SingularMatrixException extends org.ejml.data.SingularMatrixException {
		private static final long serialVersionUID = 1L;
		
		public final SimpleMatrix matrix;
		
		public SingularMatrixException(org.ejml.data.SingularMatrixException cause, SimpleMatrix matrix, String message) {
			super(message);
			this.initCause(cause);
			this.matrix = matrix;
		}
	}
	
	public static String matToRowMajorArrayString(SimpleMatrix m){
		StringBuilder sb = new StringBuilder();
		sb.append('[');
		int i;
		for(i = 0; i < m.getNumElements()-1; i++){
			double v = m.get(i);
			if(v >= 0){
				sb.append(' ');
			}
			sb.append(String.format(Locale.US,"%.3f", v));
			sb.append(',');
		}
		if(i<m.getNumElements()){
			double v = m.get(i);
			if(v >= 0){
				sb.append(' ');
			}
			sb.append(String.format(Locale.US,"%.3f", v));
		}
		sb.append(']');
		return sb.toString();
	}
	
	public static String matTo2DArrayString(SimpleMatrix m, boolean transpose){
		StringBuilder sb = new StringBuilder();
		sb.append('[');
		
		int numRows = transpose ? m.numCols():m.numRows();
		int numCols = transpose ? m.numRows():m.numCols();
		
		int r;
		for(r=0; r<numRows; r++){
			if(r!=0)
				sb.append(' ');
			int c;
			for(c=0; c<numCols-1; c++){
				double v = transpose ? m.get(c,r):m.get(r,c);
				if(v >= 0){
					sb.append(' ');
				}
				sb.append(String.format(Locale.US,"%.3f", v));
				sb.append(',');
			}
			if(c<numCols){
				double v = transpose ? m.get(c,r):m.get(r,c);
				if(v >= 0){
					sb.append(' ');
				}
				sb.append(String.format(Locale.US,"%.3f", v));
			}
			if(r != numRows-1){
				sb.append(',');
				sb.append('\n');
			}
		}
		sb.append(']');
		return sb.toString();
	}
	

	public static SimpleMatrix vectorOf(double...values){
		return vector(values.length, values);
	}
	
	public static SimpleMatrix normalize(SimpleMatrix vec, double tol){
		double norm = vec.normF();
		if(norm < tol){
			return vec;
		} else {
			return vec.scale(1.0/norm);
		}
	}
	
	public static SimpleMatrix normalize(SimpleMatrix vec){
		return normalize(vec, 0.00001);
	}
	
	
	
	public static SimpleMatrix normalizeInPlace(SimpleMatrix vec, double tol){
		double norm = vec.normF();
		if(norm < tol){
			return vec;
		} else {
			double divByNorm = 1.0/norm;
			for(int i=0; i<vec.getNumElements(); i++){
				vec.set(i, vec.get(i)*divByNorm);
			}
			return vec;
		}
	}
	
	public static SimpleMatrix normalizeInPlace(SimpleMatrix vec){
		return normalizeInPlace(vec,0.00001);
	}
	
	public static DoubleStream streamValues(SimpleMatrix vec){
		return IntStream.range(0, vec.getNumElements()).mapToDouble(vec::get);
	}
	
	public static Stream<SimpleMatrix> streamRows(SimpleMatrix mat){
		return IntStream.range(0, mat.numRows()).mapToObj(r->mat.rows(r, r+1));
	}
	
	public static Stream<SimpleMatrix> streamColumns(SimpleMatrix mat){
		return IntStream.range(0, mat.numCols()).mapToObj(c->mat.cols(c, c+1));
	}
	
	public static SimpleMatrix stackVectors(SimpleMatrix ... vecs) {
		int n = Arrays.stream(vecs).mapToInt(SimpleMatrix::getNumElements).sum();
		SimpleMatrix v = vector(n);
		int i=0;
		for(SimpleMatrix vec : vecs) {
			for(int j=0; j<vec.getNumElements(); j++) {
				v.set(i++, vec.get(j));
			}
		}
		return v;
	}
	
	public static SimpleMatrix rowmajorMat(double[][] array) {
		SimpleMatrix m = matrix(array.length, array[0].length);
		for(int i=0; i<m.numRows(); i++)
			for(int j=0; j<m.numCols(); j++) {
				m.set(i, j, array[i][j]);
			}
		return m;
	}
}
