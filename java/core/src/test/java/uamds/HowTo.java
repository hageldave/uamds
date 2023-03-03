package uamds;

import java.util.Arrays;

import org.ejml.data.DMatrixRMaj;

import hageldave.optisled.ejml.MatCalcEJML;
import hageldave.optisled.generic.numerics.MatCalc;
import uamds.other.NRV;
import uamds.other.NRVSet;
import hageldave.utils.Ref;

public class HowTo {

	public static void doLinAlg() {
		MatCalc<DMatrixRMaj> mc = new MatCalcEJML();
		DMatrixRMaj zeros = mc.zeros(3);
		DMatrixRMaj ones = mc.add(zeros, 1.0);
		DMatrixRMaj permutation = mc.matOf(new double[][] {{0,1,0},{0,0,1},{1,0,0}});
	}
	
	public static void doLinAlgGeneric() {
		MatCalc<DMatrixRMaj> mc = new MatCalcEJML();
		DMatrixRMaj zeros = mc.zeros(3);
		DMatrixRMaj ones = mc.add(zeros, 1.0);
		DMatrixRMaj permutation = mc.matOf(new double[][] {{0,1,0},{0,0,1},{1,0,0}});
	}
	
	public static <M> void doLinAlgGeneric(MatCalc<M> mc) {
		M identity = mc.eye(3);
		M scaling = mc.scale(identity, 1.0/mc.numRows(identity));
		M vec = mc.vecOf(1.0, 2.0, 3.0);
		M result = mc.mult_aTbc(vec, scaling, vec); // a^T * b * c
		double value = mc.get(result, 0, 0);
	}
	
	public static <M> void doDatasetThings(MatCalc<M> mc) {
		int numDims = 3;
		NRV<M> standardNormal = new NRV<>(mc, numDims);
		NRV<M> myNRV = new NRV<>(mc, mc.vecOf(1.0, 2.0, 3.0), mc.eye(numDims, 1.337));
		NRV<M> randCovNRV = new NRV<>(mc, mc.zeros(numDims), NRV.randCov(mc, numDims));
		
		NRVSet<M> dataset = new NRVSet<>();
		dataset.addAll(Arrays.asList(standardNormal, myNRV, randCovNRV));
		M dataSetSamples = dataset.stream()
			// draw 1000 samples form each random vector
			.map(nrv -> nrv.drawSamples(1000))
			// stack all samples to create a large matrix
			.reduce(mc::concatVert)
			.get();
		double[][] samples = mc.toArray2D(dataSetSamples);
	}
	
	public static <M> void doUAMDS(MatCalc<M> mc, NRVSet<M> dataset) {
		int lowDims = 2;
		UAMDS<M> uamds = new UAMDS<>(mc, lowDims);
		NRVSet<M> projectedDataset = uamds.calculateProjection(dataset, null, null);
		
		Ref<M[][]> affineTransf = new Ref<>();
		NRVSet<M> intermediateResult = uamds.calculateProjection(dataset, null, affineTransf);
		M[][] initialization = affineTransf.get();
		NRVSet<M> refinedResult = uamds.calculateProjection(dataset, initialization, affineTransf);
		
		int nIterations = 200;
		Ref<double[][]> stress = new Ref<>();
		Ref<double[][][]> stressDetailed = new Ref<>();
		uamds.setStochasticGDEnabled(true);
		uamds.gd.maxLineSearchIter = 10;
		uamds.calculateProjection(dataset, affineTransf.get(), affineTransf, nIterations, stress, stressDetailed);
		
		M[][] init1 = Initialization.initRandom(mc, dataset);
		M[][] init2 = Initialization.initFromUAPCA(mc, dataset);
		M[][] init3 = Initialization.initWithoutVariance(mc, dataset);
	}

}
