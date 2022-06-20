package demo;

import java.awt.Color;
import java.util.ArrayList;

import datasets.StudentGrades;
import optimization.ejml.MatCalcEJML;
import optimization.generic.numerics.MatCalc;
import uamds.NRV;
import uamds.NRVSet;
import uamds.Ref;
import uamds.UAMDS;

public class Example {

	public static void main(String[] args) {
		/* create wrapper for linear algebra lib of choice (matrix calculator) */
		MatCalcEJML mc = new MatCalcEJML();
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
	
	public static <M> void executeExample(MatCalc<M> mc) {
		/* prepare data */
		NRVSet<M> data = getData_StudentGrades(mc);
		// RVPointSet<M> data = getData_RandomizedDistribs(mc);
		
		/* prepare objects for UAMDS */
		UAMDS<M> uamds = new UAMDS<>(mc,2);
		M[][] init = null;
		Ref<M[][]> result = new Ref<>();
		Ref<double[][]> pairwiseLoss = new Ref<>();
		/* perform 10 iterations of UAMDS */
		NRVSet<M> projectedData = uamds.calculateProjection(
				data, 
				init, 
				result,
				10, // number of descend steps 
				pairwiseLoss);
		
		/* report on current loss */
		double totalLoss = 0;
		for(int i=0; i<data.size(); i++)
			for(int j=i; j<data.size(); j++)
				totalLoss += pairwiseLoss.get()[i][j];
		System.out.format("total loss : %.3f%n", totalLoss);
		System.out.println("------------------------------");
		
		/* perform another 2000 iteration of UAMDS (20 x 100 iterations) */
		for(int k=0; k<20; k++) {
			init = result.get(); // use previous result as initialization
			projectedData = uamds.calculateProjection(
					data, 
					init, 
					result,
					100, // number of descend steps 
					pairwiseLoss);

			/* report on current loss */
			totalLoss = 0;
			for(int i=0; i<data.size(); i++)
				for(int j=i; j<data.size(); j++)
					totalLoss += pairwiseLoss.get()[i][j];
			System.out.format("total loss : %.3f%n", totalLoss);
			System.out.println("------------------------------");
		}
		
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
