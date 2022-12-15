package uamds;

import java.util.Arrays;
import java.util.stream.Collectors;

import uamds.optimization.generic.numerics.MatCalc;
import uamds.other.NRV;
import uamds.other.NRVSet;
import uamds.other.Ref;
import uamds.other.UAPCA;

/**
 * Helper class to create initializations for {@link UAMDS}.
 * UAMDS has an optional initialization argument M[][] which is
 * an initial guess of the solution, i.e., the local projection matrices B_i
 * and low-dimensional means c_i.
 * <p>
 * This class offers different kinds of initializations strategies.
 */
public class Initialization {
	
	/**
	 * Performs Uncertainty-Aware PCA on the data and returns the global projection 
	 * as initialization for UAMDS.
	 * @param <M> matrix data type
	 * @param mc matrix calculator 
	 * @param data dataset
	 * @return initialization for {@link UAMDS} in matrix array form {B[],c[],P[],t[]} 
	 * as expected by {@link UAMDS#calculateProjection(NRVSet, Object[][], Ref, int, Ref, Ref)}.
	 */
	public static <M> M[][] initFromUAPCA(MatCalc<M> mc, NRVSet<M> data) {
		M[] nrvBases = data.stream().map(nrv->mc.svd(nrv.cov, true)[0]).toArray(mc::matArray);
		M[] uapca = new UAPCA<>(mc).calculate(data);
		M initProj = mc.trp(uapca[0]);
		// extract 2D projection
		initProj = mc.getRange(initProj, 0, 2, 0, mc.numCols(initProj));
		M[][] init = mc.matArray(4,data.size());
		for(int i=0; i<data.size(); i++) {
			init[0][i] = mc.mult_ab(initProj,nrvBases[i]);
			init[1][i] = mc.mult_ab(initProj,data.get(i).mean);
			init[2][i] = initProj;
			init[3][i] = mc.zeros(2);
		}
		return init;
	}
	
	/**
	 * Performs UAMDS on the data without variance and returns the resulting projections 
	 * as initialization for UAMDS.
	 * @param <M> matrix data type
	 * @param mc matrix calculator 
	 * @param data dataset
	 * @return initialization for {@link UAMDS} in matrix array form {B[],c[],P[],t[]} 
	 * as expected by {@link UAMDS#calculateProjection(NRVSet, Object[][], Ref, int, Ref, Ref)}.
	 */
	public static <M> M[][] initWithoutVariance(MatCalc<M> mc, NRVSet<M> data) {
		NRVSet<M> novarianceData = data.stream()
				.map(NRV::copy)
				.map(nrv->{
					mc.scale_inp(nrv.cov,0); 
					return nrv;
					})
				.collect(Collectors.toCollection(NRVSet::new));
		Ref<M[][]> proj = new Ref<>();
		new UAMDS<>(mc).calculateProjection(novarianceData, initRandom(mc, novarianceData), proj, 1000);
		return proj.get();
	}
	
	/**
	 * Returns randomly initialized projection matrices (B[]) and low dimensional means (c[]).
	 * The means are scaled by the median distance between dataset means, 
	 * to ensure appropriate order of magnitude of distances in low-dimensional space.
	 * @param <M> matrix data type
	 * @param mc matrix calculator 
	 * @param data dataset
	 * @return initialization for {@link UAMDS} in matrix array form {B[],c[],P[],t[]} 
	 * as expected by {@link UAMDS#calculateProjection(NRVSet, Object[][], Ref, int, Ref, Ref)}.
	 */
	public static <M> M[][] initRandom(MatCalc<M> mc, NRVSet<M> data) {
		int n = data.size();
		int loDim = 2;
		int hiDim = data.get(0).d;
		M[] nrvBases = data.stream().map(nrv->mc.svd(nrv.cov, true)[0]).toArray(mc::matArray);
		
		M pdMeans = pdMeans(mc, data);
		double[] distances = mc.toArray(pdMeans);
		Arrays.sort(distances);
		double medianDist = distances[(n+n*n)/2];
		
		M[] B = mc.matArray(n);
		M[] c = mc.matArray(n);
		for(int i=0;i<B.length; i++) {
			B[i] = mc.sub(mc.rand(loDim, hiDim),0.5);
			c[i] = mc.scale(mc.sub(mc.rand(loDim),0.5),medianDist);
		}
		/* from random B,c calculate corresponding P,t */
		M[][] init = mc.matArray(4,n);
		for(int i=0; i<data.size(); i++) {
			init[0][i] = B[i];
			init[1][i] = c[i];
			init[2][i] = mc.mult_abT(B[i], nrvBases[i]);
			init[3][i] = mc.sub(c[i], mc.mult_ab(init[2][i],data.mean(i)));
		}
		return init;
	}
	
	private static <M> M pdMeans(MatCalc<M> mc, NRVSet<M> data){
		M means = data.stream()
		.map(nrv->nrv.mean)
		.map(mc::trp)
		.reduce(mc::concatVert)
		.get();
		return mc.sqrt_inp(mc.pairwiseDistances2(means, means));
	}
	
}
