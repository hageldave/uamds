package uamds;

import java.util.Arrays;
import java.util.stream.Collectors;

import optimization.generic.numerics.MatCalc;

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
	public M[] calculate(RVPointSet<M> data) {
		/* calculating (uncertainty aware) covariance matrix for Eigen decomp */
		
		// empirical mean
		M mu = data.stream().map(nrv->nrv.mean).reduce(mc::add).get();
		mu = mc.scale_inp(mu, 1.0/data.size());
		M centering = mc.matmul(mu, mc.trp(mu));
		// average covariance
		M avgCov = data.stream().map(nrv->nrv.cov).reduce(mc::add).get();
		avgCov = mc.scale_inp(avgCov, 1.0/data.size());
		// sample covariance
		M meanMatrix = data.stream().map(nrv->mc.trp(nrv.mean)).reduce(mc::concatVert).get();
		M sampleCov = NRV.estimateFromData(mc, meanMatrix).cov;
		
		M uaCov = mc.sub( mc.add(sampleCov, avgCov), centering );
		// eigendecomposition yields principal vectors (SVD is equivalent to EVD here, but yields ordered vectors)
		M[] usv = mc.svd(uaCov, true);
		return Arrays.copyOf(usv, 2);
	}
	
	public RVPointSet<M> projectData(RVPointSet<M> data, int d) {
		return projectData(data, d, null);
	}
	
	public RVPointSet<M> projectData(RVPointSet<M> data, int d, Ref<M> proj) {
		M[] pca = calculate(data);
		M projection = mc.getRange(pca[0], 0, mc.numRows(pca[0]), 0, d);
		M transform = mc.trp(projection);
		
		RVPointSet<M> projected = data.stream()
				.map(nrv->nrv.transform(transform))
				.collect(Collectors.toCollection(RVPointSet::new));
		if(proj != null)
			proj.set(transform);
		return projected;
	}
	
}
