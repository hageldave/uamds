package uamds.other;

import java.io.Serializable;

import uamds.optimization.generic.numerics.MatCalc;

/**
 * Normally distributed random vector
 * 
 * @param <M> matrix data type
 */
public class NRV<M> implements Serializable {
	private static final long serialVersionUID = 1L;
	
	/** dimensionality */
	public int d;
	/** mean (column vector) */
	public M mean;
	/** covariance matrix (positive semi-definite) */
	public M cov;
	/** matrix calculator for generic matrix data type */
	public MatCalc<M> mc;
	
	/**
	 * standard normal
	 * @param d dimensionality
	 */
	public NRV(MatCalc<M> mc, int d) {
		this(mc, mc.zeros(d), mc.eye(d));
	}
	
	/**
	 * normally distributed random vector with specified mean and covariance
	 * @param mc
	 * @param mean
	 * @param cov
	 */
	public NRV(MatCalc<M> mc, M mean, M cov) {
		Utils.requireEquals(mc.numElem(mean)*mc.numElem(mean), mc.numElem(cov), ()->"covariance size incompatible with mean size");
		this.mc=mc;
		this.d=mc.numElem(mean);
		this.mean= mc.numCols(mean)==1 ? mean:mc.trp(mean);
		this.cov=cov;
	}
	
	/** subtract */
	public NRV<M> sub(NRV<M> other) {
		return new NRV<>(mc, mc.sub(mean, other.mean), mc.add(cov, other.cov));
	}
	
	/** mean = mean+m*delta */
	public NRV<M> translateMean(M deltaMean, double m) {
		return new NRV<>(mc, mc.add(mean, mc.scale(deltaMean, m)), cov);
	}
	
	/** given transform A, returns NRV{ A*mean, A*cov*A' } */
	public NRV<M> transform(M transform) {
		Utils.requireEquals(mc.numCols(transform), d, 
				()->"transformation does not have required number of columns so that meanNew=transform*mean is applicable. cols=" + mc.numCols(transform) + " required="+d);
		M mean = mc.matmul(transform, this.mean);
		M cov = mc.matmul(transform, mc.matmul(this.cov, mc.trp(transform)));
		return new NRV<>(mc, mean, cov);
	}
	
	/**
	 * given transform A and translation b, returns NRV{ A*mean+b, A*cov*A' }
	 * @param transform matrix
	 * @param translate vector but applied after transform
	 * @return
	 */
	public NRV<M> affineTransform(M transform, M translate) {
		NRV<M> transformed = transform(transform);
		return transformed.translateMean(translate, 1.0);
	}
	
	public NRV<M> copy() {
		return new NRV<>(mc, mc.copy(mean), mc.copy(cov));
	}
	
	public NRV<M> scaleCov(double s){
		NRV<M> copy = copy();
		mc.scale_inp(copy.cov, s);
		return copy;
	}
	
	/** probability density function */
	public MultivariateGaussian<M> pdf() {
		return MultivariateGaussian.fromNRV(this);
	}
	
	/**
	 * draws n samples from the multivariate normal distribution described by this random vector.
	 * @param n number of samples
	 * @return a matrix of n rows by d columns where each row vector is a sample.
	 */
	public M drawSamples(int n) {
		M samples = mc.randN(d, n);
		M U = mc.trp(mc.cholesky(cov));
		return mc.trp( mc.addColVec(mc.matmul(U, samples), mean) );
	}
	
	/**
	 * draws n samples from the multivariate normal distribution described by this random vector.
	 * For degenerate covariance matrices (where some singular values are zero), 
	 * the regularization parameter can be used to add a small constant diagonal to the covariance.
	 * <br>
	 * {@code cov' = (cov + I*r)}
	 * 
	 * @param n number of samples
	 * @param regularization parameter r (small,positive)
	 * @return a matrix of n rows by d columns where each row vector is a sample.
	 */
	public M drawSamples(int n, double regularization) {
		M samples = mc.randN(d, n);
		M U = mc.trp(mc.cholesky(mc.add(cov, mc.eye(d, regularization))));
		return mc.trp( mc.addColVec(mc.matmul(U, samples), mean) );
	}
	
	@Override
	public String toString() {
		return String.format("NRV{mu=%s :: sigma=%s}", Utils.arrayToString(mc.toArray(mean), 3), Utils.arrayToString(mc.toArray(cov),3));
	}
	
	/**
	 * @return string containing Java code to create this NRV
	 */
	public String toConstructorString() {
		return String.format("new NRV<>(mc, mc.vecOf(%s), mc.matOf(%d, %s))", 
				Utils.arrayToString(mc.toArray(mean), 3).replace('[', ' ').replace(']', ' ').trim(),
				d,
				Utils.arrayToString(mc.toArray(cov), 3).replace('[', ' ').replace(']', ' ').trim()
				);
	}
	
	/**
	 * @return axes of distribution (ellipse axes) as columns
	 */
	public M axes() {
		M[] svd = mc.svd(cov, true);
		return mc.mult_aTb(svd[2], mc.sqrt_inp(svd[1]));
	}

	/**
	 * creates a random covariance matrix (positive semi-definite)
	 * @param mc matrix calculator for generic matrix data type
	 * @param d dimensionality
	 * @param <M> matrix data type
	 * @return a covariance matrix
	 */
	public static <M> M randCov(MatCalc<M> mc, int d) {
		M rand = mc.sub(mc.scale(mc.rand(d, d), 2), 1);
		rand = mc.add(rand, mc.trp(rand));
		M[] evd = mc.sortEVD_inp(mc.symEvd(rand));
		double minEig = mc.get(evd[1], d*d-1);
		if(minEig < 0) {
			rand = mc.sub(rand, mc.eye(d, minEig*1.1));
		}
		if(mc.det(rand) <= 0 || mc.frob(mc.sub(mc.trp(rand), rand)) > 0.0001) {
			System.err.println("oooopsiiee");
		}
		return rand;
	}
	
	/**
	 * creates a random NRV
	 * @param mc matrix calculator for generic matrix data type
	 * @param d dimensionality
	 * @param <M> matrix data type
	 * @return
	 */
	public static <M> NRV<M> randNRV(MatCalc<M> mc, int d) {
		return new NRV<>(mc, mc.scale(mc.sub(mc.rand(d), .5), 2), mc.scale(randCov(mc,d), 1));
	}
	
	/**
	 * Estimates a multivariate normal distribution from the specified data.
	 * @param data matrix where each row vector is a data sample
	 * @return estimated normally distributed random vector
	 */
	public static <M> NRV<M> estimateFromData(MatCalc<M> mc, M data) {
		M mean = mc.colMeans(data);
		data = mc.subRowVec(data, mean);
		M cov = mc.matmul(mc.trp(data), data);
		cov = mc.scale_inp(cov, 1.0/mc.numRows(data));
		
		return new NRV<M>(mc, mc.trp(mean), cov);
	}
	
	/**
	 * Estimates a multivariate normal distribution from the specified data.
	 * @param data where each row ({@code data[rowindex][colindex]}) is a data sample
	 * @return estimated normally distributed random vector
	 */
	public static <M> NRV<M> estimateFromData(MatCalc<M> mc, double[][] data){
		return estimateFromData(mc, mc.matOf(data));
	}
	
}
