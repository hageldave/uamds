package uamds.demo;

import java.awt.Color;
import java.util.ArrayList;

import uamds.UAMDS;
import uamds.datasets.StudentGrades;
import uamds.optimization.ejml.MatCalcEJML;
import uamds.optimization.generic.numerics.MatCalc;
import uamds.other.NRV;
import uamds.other.NRVSet;
import uamds.other.Ref;

/**
 * This example demonstrates how to perform UAMDS on a data set consisting
 * of multivariate normal distributions.
 * All methods use a generic matrix data type 'M' which depends on the linear
 * algebra library used.
 * In this example the EJML library is used and all matrix operations are
 * performed via {@link MatCalcEJML} which infers the type DMatrixRMaj for 'M'.
 */
public class Example {

	public static void main(String[] args) {
		/* create wrapper for linear algebra lib of choice (matrix calculator) */
		MatCalc<?> mc = new MatCalcEJML();
		/* run example */
		executeExample(mc);
	}
	
	public static <M> NRVSet<M> getData_StudentGrades(MatCalc<M> mc){
		return NRVSet.fromArray(StudentGrades.get(mc, 0));
	}
	
	public static <M> NRVSet<M> getData_RandomizedDistribs(MatCalc<M> mc){
		NRVSet<M> data = new NRVSet<>();
		int dimensionality = 4;
		int numInstances = 5;
		/* we are generating a number of random Gaussians for demonstration */
		for(int i=0; i<numInstances; i++) {
			M mean = mc.scale(mc.rand(dimensionality),10);
			M cov = NRV.randCov(mc, dimensionality);
			data.add(new NRV<M>(mc, mean, cov));
		}
		return data;
	}
	
	public static <M> NRVSet<M> performUAMDS(
			MatCalc<M> mc, 
			NRVSet<M> data, 
			M[][] init, 
			Ref<M[][]> result, 
			int numIters, 
			Ref<double[][]> pairwiseLoss,
			boolean stochasticGD) 
	{
		UAMDS<M> uamds = new UAMDS<>(mc,2);
		uamds.setStochasticGDEnabled(stochasticGD);
		/* perform a number of UAMDS iterations */
		NRVSet<M> projectedData = uamds.calculateProjection(
				data, 
				init, 
				result,
				numIters,
				pairwiseLoss);
		/* report on current loss */
		double totalLoss = 0;
		for(int i=0; i<data.size(); i++)
			for(int j=i; j<data.size(); j++)
				totalLoss += pairwiseLoss.get()[i][j];
		System.out.format("total loss : %.3f%n", totalLoss);
		System.out.println("------------------------------");
		
		return projectedData;
	}

	public static <M> void executeExample(MatCalc<M> mc) {
		/* prepare data */
		NRVSet<M> data = getData_StudentGrades(mc);
		// RVPointSet<M> data = getData_RandomizedDistribs(mc);
		
		/* prepare objects for UAMDS */
		M[][] init = null;
		Ref<M[][]> result = new Ref<>();
		Ref<double[][]> pairwiseLoss = new Ref<>();
		/* perform 10 iterations of UAMDS */
		performUAMDS(mc, data, init, result, 10, pairwiseLoss, false);
		/* perform another 2000 iteration of UAMDS (20 x 100 iterations) with stochastic gradient descent */
		for(int k=0; k<20; k++) {
			init = result.get(); // use previous result as initialization
			performUAMDS(mc, data, init, result, 100, pairwiseLoss, true);
		}
		/* perform some final iterations with regular gradient descent */
		init = result.get(); // use previous result as initialization
		NRVSet<M> projectedData = performUAMDS(mc, data, init, result, 100, pairwiseLoss, false);
		
		/* report on final pairwise loss */
		for(int i=0; i<data.size(); i++)
			for(int j=i; j<data.size(); j++)
				System.out.format("loss between %d and %d : %.3f%n", i,j,pairwiseLoss.get()[i][j]);
		System.out.println("------------------------------");
		/* print projected random vectors */
		for(NRV<M> nrv : projectedData) {
			System.out.println(nrv);
		}
		/* low fidelity visualization */
		LoFiScatter scatter = new LoFiScatter();
		scatter.setBackground(Color.white);
		scatter.display("UAMDS demo");
		/* draw samples from projected distributions */
		ArrayList<double[][]> sampleSet = new ArrayList<>();
		for(NRV<M> nrv : projectedData) {
			M samples = nrv.drawSamples(1000, 0.01);
			sampleSet.add(mc.toArray2D(samples));
		}
		scatter.setPointSets(sampleSet);
		scatter.setPointOpacity(255/3);
	}
	
}
