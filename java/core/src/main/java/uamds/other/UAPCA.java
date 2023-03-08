package uamds.other;

import java.util.Arrays;
import java.util.stream.Collectors;

import hageldave.optisled.generic.numerics.MatCalc;
import hageldave.utils.Ref;

/**
 * Implementation of "uncertainty-aware principal component analysis"
 * @see <a href="https://doi.org/10.1109/TVCG.2019.2934812">Publication</a>
 * @param <M> matrix type
 */
public class UAPCA<M> {
	
	final MatCalc<M> mc;
	
	public UAPCA(MatCalc<M> mc) {
		this.mc = mc;
	}

	/**
	 * Calculates uncertainty aware PCA for a set of random vectors
	 * @param data set of random vectors
	 * @return projection matrix (each column is a principal vector), and vector of eigenvalues
	 */
	public M[] calculate(NRVSet<M> data) {
		M uaCov = calcUACov(data);
		// eigendecomposition yields principal vectors (SVD is equivalent to EVD here, but yields ordered vectors)
		M[] usv = mc.svd(uaCov, true);
		return Arrays.copyOf(usv, 2);
	}
	
	public M calcUACov(NRVSet<M> data) {
		/* calculating (uncertainty-aware) covariance matrix for Eigen decomp */
		
		// empirical mean
		M mu = data.stream().map(nrv->nrv.mean).reduce(mc::add).get();
		mu = mc.scale_inp(mu, 1.0/data.size());
		M centering = mc.matmul(mu, mc.trp(mu));
		// average covariance
		M avgCov = data.stream().map(nrv->nrv.cov).reduce(mc::add).get();
		avgCov = mc.scale_inp(avgCov, 1.0/data.size());
		// sample covariance
		M sampleCov = data.stream().map(nrv->mc.mult_abT(nrv.mean, nrv.mean)).reduce(mc::add).get();
		mc.scale_inp(sampleCov, 1.0/data.size());
		
		M uaCov = mc.sub( mc.add(sampleCov, avgCov), centering );
		return uaCov;
	}
	
	public M calcUACovV2(NRVSet<M> data) {
		/* calculating (uncertainty-aware) covariance matrix for Eigen decomp */
		final double n = data.size();
		// empirical mean
		M mu = data.stream().map(nrv->nrv.mean).reduce(mc::add).get();
		mu = mc.scale_inp(mu, 1.0/n);
		M centering = mc.mult_abT(mu, mu);
		// average covariance
		M avgCov = data.stream().map(nrv->nrv.cov).reduce(mc::add).get();
		avgCov = mc.scale_inp(avgCov, 1.0/n);
		// sample covariance
		M sampleCov = data.stream().map(nrv->mc.mult_abT(nrv.mean, nrv.mean)).reduce(mc::add).get();
		mc.scale_inp(sampleCov, 1.0/n);
		// inter-distribution covariance
		M interCov = data.stream().map(nrv->nrv.cov).reduce(mc::add).get();
		interCov = mc.add(interCov, 0); // placeholder, here goes the sum of covariances between distributions
		mc.scale_inp(interCov, 1.0/(n*n));
		
		M pos = mc.add(sampleCov, avgCov);
		M neg = mc.add(centering, interCov);
		M uaCov = mc.sub(pos, neg);
		return uaCov;
	}
	
	public NRVSet<M> projectData(NRVSet<M> data, int d) {
		return projectData(data, d, null);
	}
	
	public NRVSet<M> projectData(NRVSet<M> data, int d, Ref<M> proj) {
		M[] pca = calculate(data);
		M projection = mc.getRange(pca[0], 0, mc.numRows(pca[0]), 0, d);
		M transform = mc.trp(projection);
		
		NRVSet<M> projected = data.stream()
				.map(nrv->nrv.transform(transform))
				.collect(Collectors.toCollection(NRVSet::new));
		if(proj != null)
			proj.set(transform);
		return projected;
	}
	
}
