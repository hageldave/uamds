package uamds.other;

import uamds.optimization.generic.numerics.MatCalc;

/**
 * Gaussian probability density function
 *
 * @param <M> matrix data type
 */
public class MultivariateGaussian<M> {

	/** inverse of covariance matrix */
	public final M covInv;
	/** mean */
	public final M mu;
	/** scaling term: 1 / sqrt((2π)^d |Σ|) */
	public final double scaling;
	/** number of dimensions */
	public final int d;
	/** matrix calculator for generic matrix data type */
	public final MatCalc<M> mc;
	
	
	public MultivariateGaussian(MatCalc<M> mc, M covInv, M mu, double scaling) {
		Utils.requireEquals(mc.numRows(covInv), mc.numCols(covInv), ()->"nonsquare matrix");
		Utils.requireEquals(mc.numCols(covInv), mc.numElem(mu), ()->"incompatible dimensions between mu and covInv");
		this.mc = mc;
		this.covInv = covInv;
		this.mu = mc.numCols(mu)==1 ? mu:mc.trp(mu);
		this.scaling = scaling;
		this.d = mc.numElem(mu);
	}
	
	public static <M> MultivariateGaussian<M> fromNRV(NRV<M> nrv) {
		MatCalc<M> mc = nrv.mc;
		M cov = mc.copy(nrv.cov);
		double det = mc.det(cov);
		if(det<0) System.err.println("negative determinant: " + det + " for " + nrv);
		M[] svd = mc.svd(cov, true);
		for(int i=0;i<nrv.d;i++)
			mc.set_inp(svd[1], i, i, 1.0/mc.get(svd[1], i, i)); // inverting S
		M covInv = mc.mult_ab(svd[2], mc.mult_abT(svd[1], svd[2]));
		double scaling = (Math.pow(2*Math.PI, -0.5*nrv.d) * Math.pow(det, -0.5));
		return new MultivariateGaussian<>(mc, covInv, mc.copy(nrv.mean), scaling);
	}
	
	public double evalAt(M x) {
		var diff = mc.sub(x,mu);
		double exponent = -0.5 *  mc.get( mc.matmul( mc.trp(diff), mc.matmul(covInv, diff) ) , 0);
		return scaling * Math.exp(exponent);
	}
	
	public double evalAt(double...coords) {
		return evalAt(mc.vecOf(coords));
	}
	
	/* points are row vectors */
	public M evalManyAt(M points) {
		M diff = mc.subRowVec(points, mu);
		M left = mc.matmul(diff, covInv);
		M exponents = mc.elmmul(left, diff);
		exponents = mc.rowSums(exponents);
		exponents = mc.scale_inp(exponents, -0.5);
		
		return mc.scale_inp(mc.exp_inp(exponents), scaling);
	}
	
	
}
