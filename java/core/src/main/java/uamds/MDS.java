//package uamds;
//
//import static uamds.Utils.sq;
//
//import java.util.ArrayList;
//import java.util.List;
//import java.util.stream.Collectors;
//import java.util.stream.IntStream;
//
//import org.ejml.data.DMatrixRMaj;
//
//import hageldave.optisled.ejml.MatCalcEJML;
//import hageldave.optisled.generic.numerics.MatCalc;
//import hageldave.optisled.generic.numerics.NumericGradient;
//import hageldave.optisled.generic.problem.ScalarFN;
//import hageldave.optisled.generic.problem.VectorFN;
//import hageldave.optisled.generic.solver.GradientDescent;
//
//public class MDS<M> {
//	
//	final MatCalc<M> mc;
//	
//	public MDS(MatCalc<M> mc) {
//		this.mc = mc;
//	}
//
//	public RVPointSet<M> monteCarlo(RVPointSet<M> data, int samplesPerDistribution) {
//		List<M> samples = data.stream().map(nrv->nrv.drawSamples(samplesPerDistribution)).collect(Collectors.toList());
//		M allSamples = samples.stream().reduce(mc::concatVert).get();
//		int maxDescendSteps = mc.numRows(allSamples)*5;
//		M low = this.calculate(allSamples, null, maxDescendSteps);
//		RVPointSet<M> lowDistributions = new RVPointSet<>();
//		for(int i=0; i<data.size(); i++) {
//			M lowSamples = mc.getRange(low, i*samplesPerDistribution, (i+1)*samplesPerDistribution, 0, mc.numCols(low));
//			NRV<M> nrv = NRV.estimateFromData(mc, lowSamples);
//			lowDistributions.add(nrv);
//		}
//		return lowDistributions;
//	}
//	
//	public RVPointSet<M> independentMonteCarlo(RVPointSet<M> data, RVPointSet<M> init, int samplesPerDistribution) {
//		RVPointSet<M> lowDistributions = init;
//		if(lowDistributions == null) {
//			lowDistributions = new RVPointSet<>();
//			// init low dim distributions randomly
//			for(int i=0; i<data.size(); i++) {
//				lowDistributions.add(new NRV<>(mc,mc.rand(2), NRV.randCov(mc,2)));
//			}
//		}
//		// repeat this until convergence
//		for(int step=0; step < 1; step++) {
//			// sample distributions
//			List<M> highSamples = data.stream().map(nrv->nrv.drawSamples(samplesPerDistribution)).collect(Collectors.toList());
//			List<M> lowSamples = lowDistributions.stream().map(nrv->nrv.drawSamples(samplesPerDistribution)).collect(Collectors.toList());
//			// concatenate samples into single matrix
//			M allHighSamples = highSamples.stream().reduce(mc::concatVert).get();
//			M allLowSamples = lowSamples.stream().reduce(mc::concatVert).get();
//			// do a single gradient descent step (MDS step)
//			allLowSamples = this.calculate(allHighSamples, allLowSamples, 1);
//			// estimate new low dim distributions after MDS step
//			for(int i=0; i<data.size(); i++) {
//				// extract part of matrix that belongs to samples of distribution i
//				M lowSamplesOfDistribution = mc.getRange(allLowSamples, i*samplesPerDistribution, (i+1)*samplesPerDistribution, 0, 2);
//				// estimate new distribution
//				NRV<M> nrv = NRV.estimateFromData(mc,lowSamplesOfDistribution);
//				lowDistributions.set(i,nrv);
//			}
//		}
//		return lowDistributions;
//	}
//	
//	public RVPointSet<M> independentMonteCarloCoupled(RVPointSet<M> data, RVPointSet<M> init, int samplesPerDistribution) {
//		RVPointSet<M> lowDistributions = init;
//		if(lowDistributions == null) {
//			lowDistributions = new RVPointSet<>();
//			// init low dim distributions randomly
//			for(int i=0; i<data.size(); i++) {
//				lowDistributions.add(new NRV<>(mc, mc.rand(2), NRV.randCov(mc,2)));
//			}
//		}
//		// repeat this until convergence
//		for(int step=0; step < 1; step++) {
//			// sample distributions
//			List<M> highSamples = data.stream().map(nrv->nrv.drawSamples(samplesPerDistribution)).collect(Collectors.toList());
//			List<M[]> highSVD = data.stream().map(nrv->mc.svd(nrv.cov, true)).collect(Collectors.toList());
//			List<M[]> lowSVD = lowDistributions.stream().map(nrv->mc.svd(nrv.cov, true)).collect(Collectors.toList());
//			// calculate corresponding samples in low by projecting onto principal axes
//			List<M> lowSamples = new ArrayList<>();
//			for(int i=0; i<data.size(); i++) {
//				M proj = mc.matmul(mc.subRowVec(highSamples.get(i), data.mean(i)), mc.trp(highSVD.get(i)[0]));
//				M invscale = mc.elemwise_inp(mc.diagV(highSVD.get(i)[1]), v->Math.pow(Math.sqrt(v), -1.0));
//				proj = mc.matmul(proj, mc.diagM(invscale));
//				proj = mc.getRange(proj, 0, mc.numRows(proj), 0, 2);
//				M scale = mc.sqrt_inp(mc.diagV(lowSVD.get(i)[1]));
//				proj = mc.matmul(proj, mc.diagM(scale));
//				M lowSample = mc.addRowVec(mc.matmul(proj, lowSVD.get(i)[0]), lowDistributions.mean(i));
//				lowSamples.add(lowSample);
//			}
////			DoubleMatrix[] lowSamples = lowDistributions.stream().map(nrv->nrv.drawSamples(samplesPerDistribution)).toArray(DoubleMatrix[]::new);
//			// concatenate samples into single matrix
//			M allHighSamples = highSamples.stream().reduce(mc::concatVert).get();
//			M allLowSamples = lowSamples.stream().reduce(mc::concatVert).get();
//			// do a single gradient descent step (MDS step)
//			allLowSamples = calculate(allHighSamples, allLowSamples, 1);
//			// estimate new low dim distributions after MDS step
//			for(int i=0; i<data.size(); i++) {
//				// extract part of matrix that belongs to samples of distribution i
//				M lowSamplesOfDistribution = mc.getRange(allLowSamples, i*samplesPerDistribution, (i+1)*samplesPerDistribution, 0, 2);
//				// estimate new distribution
//				NRV<M> nrv = NRV.estimateFromData(mc,lowSamplesOfDistribution);
//				lowDistributions.set(i,nrv);
//			}
//		}
//		return lowDistributions;
//	}
//
//	/**
//	 * performs MDS on the specified data (set of row vectors)
//	 * @param data highdimensional dataset, each row is a data point
//	 * @param init low dimensional initialization, each row is a data point
//	 * @return low dimensional embedding
//	 */
//	public M calculate(M data, M init, int maxDescendSteps) {
//		int n = mc.numRows(data);
//
//		M low;
//		if(init != null) {
//			low=init;
//		} else {
//			low=mc.sub(mc.scale(mc.rand(n, 2),2),1);
//		}
//		
//		/* define loss+gradient and projection function for minimizing with projected gradient descent */
//		ScalarFN<M> fx;
//		VectorFN<M> dfx;
//		
//		M x = vectorize(low, null);
//		M dx = mc.copy(x);
//		M dxAsPointSet = mc.copy(low);
//		
//		// loss function, gradient and projection mapping
//		ScalarFN<M> loss = (x_) -> stress(data, fromVectorizedForm(x_, n, low));
//		VectorFN<M> gradient = (x_) -> vectorize(stressGradient(data, fromVectorizedForm(x_, n, low), dxAsPointSet), dx);
//		NumericGradient<M> gradient_numeric = new NumericGradient<>(mc, loss); gradient_numeric.h=0.00001;
//		fx = loss;
//		dfx = (x_)->{ // analytic gradient with on the fly numeric gradient check
//			var dxCheck = gradient_numeric.central.evaluate(x_);
//			var dx_ = gradient.evaluate(x_);
//			double diff = mc.norm(mc.sub(dx_,dxCheck));
//			double lenratio = mc.norm(dx_)/mc.norm(dxCheck);
//			double dot = mc.inner(dx_,dxCheck)/(mc.norm(dx_)*mc.norm(dxCheck));
//			if(diff > 0.01)
//				System.err.println("gradient mismatch! diff="+diff+" dot="+dot+ " lenRatio="+lenratio);
//			return dx_;
//		};
////		dfx = gradient_numeric;
//		dfx = gradient;
//		
//		/* perform projected gradient descent */
//		GradientDescent<M> gradientDescent = new GradientDescent<>(mc);
//		gradientDescent.maxDescentSteps = maxDescendSteps;
//		gradientDescent.lineSearchFactor = 0.001;
//		M arg_min = gradientDescent.arg_min(fx, dfx, x, null);
//		fromVectorizedForm(arg_min, n, low);
//		System.out.println("loss="+gradientDescent.lossOnTermination);
//		
//		// move center of mass to origin
//		M center = mc.colSums(low);
//		center = mc.scale_inp(center, 1.0/mc.numRows(data));
//		center = mc.subRowVec(low, center);
//		
//		return low;
//	}
//	
//	protected M vectorize(M low, M preAlloc) {
//		if(preAlloc==null)
//			preAlloc = mc.vecOf(mc.toArray(low).clone());
//		else
//			mc.copyValues(low, preAlloc);
//		return preAlloc;
//	}
//	
//	protected M fromVectorizedForm(M v, int n, M preAlloc) {
//		if(preAlloc==null)
//			preAlloc = mc.matOf(n, mc.toArray(v).clone());
//		else
//			mc.copyValues(v, preAlloc);
//		return preAlloc;
//	}
//	
//	protected double stress(M data, M low) {
//		double sum=0;
//		for(int i=0; i<mc.numRows(data); i++) {
//			M pi = mc.getRow(data, i);
//			M xi = mc.getRow(low, i);
//			for(int j=i+1; j<mc.numRows(data); j++) {
//				M pj = mc.getRow(data, j);
//				M xj = mc.getRow(low, j);
//				M pdiff = mc.sub(pi,pj);
//				M xdiff = mc.sub(xi,xj);
//				double diffexp = mc.norm2(pdiff) - mc.norm2(xdiff);
//				sum += sq(diffexp);
//			}
//		}
//		return sum/mc.numRows(data);
//	}
//	
//	M stressGradient(M data, M low, M result) {
//		if(result==null)
//			result = mc.copy(low);
//		mc.scale(result,0);
//		
//		for(int i=0; i<mc.numRows(data); i++) {
//			M pi = mc.getRow(data, i);
//			M xi = mc.getRow(low, i);
//			for(int j=i+1; j<mc.numRows(data); j++) {
//				M pj = mc.getRow(data, j);
//				M xj = mc.getRow(low, j);
//				
//				M pdiff = mc.sub(pi, pj);
//				M xdiff = mc.sub(xi, xj);
//				
//				double diffexp = mc.norm2(pdiff) - mc.norm2(xdiff);
//				mc.scale_inp(xdiff,-diffexp);
//				for(int k=0; k<mc.numCols(low); k++) {
//					mc.set_inp(result, i, k, mc.get(result, i, k)+mc.get(xdiff, k));
//					mc.set_inp(result, j, k, mc.get(result, j, k)-mc.get(xdiff, k));
//				}
//			}
//		}
//		return mc.scale_inp(result, 4.0/mc.numRows(data));
//	}
//	
//	public static void main(String[] args) {
//		MatCalcEJML mc = new MatCalcEJML();
//		RVPointSet<DMatrixRMaj> pointset = IntStream.range(0, 5).mapToObj(i->NRV.randNRV(mc, 5)).collect(Collectors.toCollection(RVPointSet::new));
//		RVPointSet<DMatrixRMaj> mds = new MDS<>(mc).monteCarlo(pointset, 128);
//		System.out.println(mds);
//	}
//	
//}
