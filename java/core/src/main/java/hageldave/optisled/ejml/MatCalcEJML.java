package hageldave.optisled.ejml;

import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.CholeskyDecomposition_F64;
import org.ejml.interfaces.decomposition.EigenDecomposition_F64;
import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;
import org.ejml.simple.ops.SimpleOperations_DDRM;

import hageldave.optisled.generic.numerics.MatCalc;

/**
 * {@link MatCalc} implementation using {@link DMatrixRMaj} matrix type of
 * Efficient Java Matrix Library (EJML).
 */
public class MatCalcEJML implements MatCalc<DMatrixRMaj> {
	
	static final SimpleOperations_DDRM ops = new SimpleOperations_DDRM();

	@Override
	public DMatrixRMaj vecOf(double... v) {
		return DMatrixRMaj.wrap(v.length, 1, v);
	}
	
	@Override
	public DMatrixRMaj matOf(double[][] values) {
		return new DMatrixRMaj(values);
	}
	
	@Override
	public DMatrixRMaj matOf(int nRows, double... values) {
		return DMatrixRMaj.wrap(nRows, values.length/nRows, values);
	}

	@Override
	public DMatrixRMaj zeros(int size) {
		return vecOf(new double[size]);
	}
	
	@Override
	public DMatrixRMaj zeros(int rows, int columns) {
		return matOf(rows, new double[rows*columns]);
	}
	
	@Override
	public DMatrixRMaj eye(int n, double s) {
		DMatrixRMaj eye = eye(n);
		ops.scale(eye, s, eye);
		return eye;
	}
	
	@Override
	public DMatrixRMaj eye(int n) {
		DMatrixRMaj m = zeros(n, n);
		ops.setIdentity(m);
		return m;
	}

	@Override
	public int numRows(DMatrixRMaj m) {
		return m.numRows;
	}

	@Override
	public int numCols(DMatrixRMaj m) {
		return m.numCols;
	}

	@Override
	public double inner(DMatrixRMaj a, DMatrixRMaj b) {
		return ops.dot(a, b);
	}

	@Override
	public DMatrixRMaj scale(DMatrixRMaj m, double s) {
		DMatrixRMaj copy = m.copy();
		ops.scale(m, s, copy);
		return copy;
	}

	@Override
	public DMatrixRMaj scale_inp(DMatrixRMaj m, double s) {
		ops.scale(m, s, m);
		return m;
	}

	@Override
	public DMatrixRMaj matmul(DMatrixRMaj a, DMatrixRMaj b) {
		DMatrixRMaj c = zeros(a.numRows, b.numCols);
		ops.mult(a, b, c);
		return c;
	}
	
	@Override
	public DMatrixRMaj elmmul(DMatrixRMaj a, DMatrixRMaj b) {
		DMatrixRMaj c = a.copy();
		ops.elementMult(a, b, c);
		return c;
	}
	
	@Override
	public DMatrixRMaj rowSums(DMatrixRMaj m) {
		return CommonOps_DDRM.sumRows(m,null);
	}
	
	@Override
	public DMatrixRMaj colSums(DMatrixRMaj m) {
		return CommonOps_DDRM.sumCols(m,null);
	}
	
	@Override
	public DMatrixRMaj rowMins(DMatrixRMaj m) {
		return CommonOps_DDRM.minRows(m, null);
	}
	
	@Override
	public DMatrixRMaj colMins(DMatrixRMaj m) {
		return CommonOps_DDRM.minCols(m, null);
	}
	
	@Override
	public DMatrixRMaj rowMaxs(DMatrixRMaj m) {
		return CommonOps_DDRM.maxRows(m, null);
	}
	
	@Override
	public DMatrixRMaj colMaxs(DMatrixRMaj m) {
		return CommonOps_DDRM.maxCols(m, null);
	}
	
	@Override
	public DMatrixRMaj mulColsByRowVec(DMatrixRMaj m, DMatrixRMaj rowV) {
		DMatrixRMaj copy = m.copy();
		CommonOps_DDRM.multCols(copy, toArray(rowV));
		return copy;
	}
	
	@Override
	public DMatrixRMaj mulRowsByColVec(DMatrixRMaj m, DMatrixRMaj colV) {
		DMatrixRMaj copy = m.copy();
		CommonOps_DDRM.multRows(toArray(colV), copy);
		return copy;
	}
	
	@Override
	public DMatrixRMaj subRowVec(DMatrixRMaj m, DMatrixRMaj rowV) {
		int cols = m.numCols;
		int rows = m.numRows;
		DMatrixRMaj copy = m.copy();
		for(int r=0; r<rows; r++) {
			for(int c=0; c<cols; c++) {
				copy.data[r*cols+c] -= rowV.data[c];
			}
		}
		return copy;
	}
	
	@Override
	public DMatrixRMaj subColVec(DMatrixRMaj m, DMatrixRMaj colV) {
		int cols = m.numCols;
		int rows = m.numRows;
		DMatrixRMaj copy = m.copy();
		for(int r=0; r<rows; r++) {
			for(int c=0; c<cols; c++) {
				copy.data[r*cols+c] -= colV.data[r];
			}
		}
		return copy;
	}
	
	@Override
	public DMatrixRMaj addRowVec(DMatrixRMaj m, DMatrixRMaj rowV) {
		int cols = m.numCols;
		int rows = m.numRows;
		DMatrixRMaj copy = m.copy();
		for(int r=0; r<rows; r++) {
			for(int c=0; c<cols; c++) {
				copy.data[r*cols+c] += rowV.data[c];
			}
		}
		return copy;
	}
	
	@Override
	public DMatrixRMaj addColVec(DMatrixRMaj m, DMatrixRMaj colV) {
		int cols = m.numCols;
		int rows = m.numRows;
		DMatrixRMaj copy = m.copy();
		for(int r=0; r<rows; r++) {
			for(int c=0; c<cols; c++) {
				copy.data[r*cols+c] += colV.data[r];
			}
		}
		return copy;
	}
	
	@Override
	public DMatrixRMaj trp(DMatrixRMaj m) {
		DMatrixRMaj trp = zeros(m.numCols, m.numRows);
		ops.transpose(m, trp);
		return trp;
	}

	@Override
	public DMatrixRMaj add(DMatrixRMaj a, DMatrixRMaj b) {
		DMatrixRMaj c = a.copy();
		ops.plus(a, b, c);
		return c;
	}
	
	@Override
	public DMatrixRMaj add_inp(DMatrixRMaj a, DMatrixRMaj b) {
		ops.plus(a, b, a);
		return a;
	}
	
	@Override
	public DMatrixRMaj sub_inp(DMatrixRMaj a, DMatrixRMaj b) {
		ops.minus(a, b, a);
		return a;
	}
	
	@Override
	public DMatrixRMaj add(DMatrixRMaj a, double b) {
		DMatrixRMaj c = a.copy();
		ops.plus(a, b, c);
		return c;
	}

	@Override
	public DMatrixRMaj sub(DMatrixRMaj a, DMatrixRMaj b) {
		DMatrixRMaj c = a.copy();
		ops.minus(a, b, c);
		return c;
	}

	@Override
	public DMatrixRMaj set_inp(DMatrixRMaj m, int idx, double v) {
		m.set(idx, v);
		return m;
	}

	@Override
	public double get(DMatrixRMaj m, int idx) {
		return m.get(idx);
	}

	@Override
	public DMatrixRMaj copy(DMatrixRMaj m) {
		return m.copy();
	}

	@Override
	public double[] toArray(DMatrixRMaj m) {
		return m.data;
	}

	@Override
	public double[][] toArray2D(DMatrixRMaj m) {
		double[][] data = new double[m.numRows][m.numCols];
		for(int r=0; r<m.numRows; r++)
			for(int c=0; c<m.numCols; c++)
				data[r][c] = m.unsafe_get(r, c);
		return data;
	}
	
	@Override
	public DMatrixRMaj[] svd(DMatrixRMaj m, boolean full) {
		SimpleSVD<SimpleMatrix> svd = SimpleMatrix.wrap(m).svd(!full);
		return new DMatrixRMaj[] {svd.getU().getDDRM(), svd.getW().getDDRM(), svd.getV().getDDRM()};
	}
	
	@Override
	public DMatrixRMaj cholesky(DMatrixRMaj m) {
		CholeskyDecomposition_F64<DMatrixRMaj> chol = DecompositionFactory_DDRM.chol(false);
		chol.decompose(m);
		return chol.getT(null);
	}
	
	@Override
	public DMatrixRMaj[] symEvd(DMatrixRMaj m) {
		m = m.copy();
		EigenDecomposition_F64<DMatrixRMaj> eig = DecompositionFactory_DDRM.eig(true, true);
		eig.decompose(m);
		double[] vals = new double[m.numCols];
		for(int i=0; i<vals.length; i++) {
			vals[i] = eig.getEigenvalue(i).real;
		}
		DMatrixRMaj vectors = zeros(m.numCols, m.numCols);
		
		for(int i=0; i<m.numCols; i++) {
			DMatrixRMaj ev = eig.getEigenVector(i);
			ops.setColumn(vectors, i, 0, ev.data);
		}
		return new DMatrixRMaj[] {vectors, diagM(vecOf(vals))};
	}

	@Override
	public double det(DMatrixRMaj m) {
		return ops.determinant(m);
	}
	
	@Override
	public DMatrixRMaj exp_inp(DMatrixRMaj m) {
		ops.elementExp(m, m);
		return m;
	}
	
	@Override
	public DMatrixRMaj concatHorz(DMatrixRMaj a, DMatrixRMaj b) {
		return SimpleMatrix.wrap(a).concatColumns(SimpleMatrix.wrap(b)).getDDRM();
	}
	
	@Override
	public DMatrixRMaj concatVert(DMatrixRMaj a, DMatrixRMaj b) {
		return SimpleMatrix.wrap(a).concatRows(SimpleMatrix.wrap(b)).getDDRM();
	}
	
	@Override
	public DMatrixRMaj getRange(DMatrixRMaj m, int ra, int rb, int ca, int cb) {
		DMatrixRMaj dst = zeros(rb-ra, cb-ca);
		CommonOps_DDRM.extract(m, ra, rb, ca, cb, dst, 0, 0);
		return dst;
	}
	
	@Override
	public DMatrixRMaj diagM(DMatrixRMaj v) {
		if(!(v.numCols==1 || v.numRows==1)) {
			throw new IllegalArgumentException("argument has to be a vector");
		}
		return CommonOps_DDRM.diag(v.data);
	}
	
	@Override
	public DMatrixRMaj diagV(DMatrixRMaj m) {
		if(m.numCols != m.numRows) {
			throw new IllegalArgumentException("argument has to be a square matrix");
		}
		DMatrixRMaj v = zeros(m.numCols);
		CommonOps_DDRM.extractDiag(m, v);
		return v;
	}
	
	@Override
	public DMatrixRMaj[] matArray(int n) {
		return new DMatrixRMaj[n];
	}
	
	@Override
	public DMatrixRMaj[][] matArray(int m, int n) {
		return new DMatrixRMaj[m][n];
	}
	
	@Override
	public DMatrixRMaj mult_aTb(DMatrixRMaj a, DMatrixRMaj b) {
		DMatrixRMaj c = zeros(a.numCols, b.numCols);
		ops.multTransA(a, b, c);
		return c;
	}
	
	@Override
	public DMatrixRMaj mult_abT(DMatrixRMaj a, DMatrixRMaj b) {
		DMatrixRMaj c = zeros(a.numRows, b.numRows);
		CommonOps_DDRM.multTransB(a, b, c);
		return c;
	}
	
	@Override
	public void copyValues(DMatrixRMaj src, DMatrixRMaj target) {
		System.arraycopy(src.data, 0, target.data, 0, Math.min(src.data.length, target.data.length));
	}
	
	@Override
	public void copyValues(DMatrixRMaj src, int startSrc, DMatrixRMaj target, int startTarget, int len) {
		System.arraycopy(src.data, startSrc, target.data, startTarget, len);
	}
	
	@Override
	public double sum(DMatrixRMaj m) {
		double sum=0;
		for(int i=0; i<m.data.length; i++) {
			sum += m.data[i];
		}
		return sum;
	}
	
	@Override
	public double frob2(DMatrixRMaj a) {
		double total = 0;
        double scale = CommonOps_DDRM.elementMaxAbs(a);
        if (scale == 0.0)
            return 0.0;
        
        final int size = a.getNumElements();
        double divByScale = 1.0/scale;
        for (int i = 0; i < size; i++) {
            double val = a.data[i]*divByScale;
            total += val*val;
        }
        return scale*scale*total;
	}

	@Override
	public DMatrixRMaj pinv(DMatrixRMaj m) {
		DMatrixRMaj result = zeros(numCols(m), numRows(m));
		ops.pseudoInverse(m, result);
		return result;
	}
}
