package experiments;

import java.awt.Color;
import java.util.ArrayList;

import uamds.datasets.StudentGrades;
import uamds.demo.LoFiScatter;
import uamds.optimization.ejml.MatCalcEJML;
import uamds.optimization.generic.numerics.MatCalc;
import uamds.other.NRV;
import uamds.other.NRVSet;
import uamds.other.UAPCA;

public class UAPCATest {

	public static void main(String[] args) {
		test(new MatCalcEJML());
	}
	
	static <M> void test(MatCalc<M> mc) {
//		NRVSet<M> dataset = new NRVSet<>();
//		
//		dataset.add( new NRV<>(mc, 
//				mc.vecOf(0,0), mc.diagM(mc.vecOf(2,1)))
//				);
//		dataset.add( new NRV<>(mc, 
//				mc.vecOf(18,18), mc.diagM(mc.vecOf(16,2)))
//				);
////		dataset.add( new NRV<>(mc, 
////				mc.vecOf(-2,7), mc.diagM(mc.vecOf(1,12)))
////				);
//		dataset.add( new NRV<>(mc, 
//				mc.vecOf(-2,7), cov2Drot(mc, mc.vecOf(1.0, 8.0), Math.PI/4))
//				);
		
		NRVSet<M> dataset = NRVSet.fromArray(StudentGrades.get(mc, 0));
		{
			// sanitize covariance matrices
			dataset.forEach(nrv->mc.add_inp(nrv.cov, mc.eye(mc.numCols(nrv.cov), 0.0001)));
			// center means
			M meanofmeans = mc.scale(dataset.stream().map(nrv->nrv.mean).reduce(mc::add).get(), 1.0/dataset.size());
			dataset.forEach(nrv->mc.sub_inp(nrv.mean, meanofmeans));
		}
		
		final int d = dataset.get(0).d;
		
		LoFiScatter scatter = new LoFiScatter();
		ArrayList<double[][]> sampleSet = new ArrayList<>();
		for(NRV<M> nrv : dataset) {
			M samples = nrv.drawSamples(1000);
			sampleSet.add(mc.toArray2D(samples));
		}
		scatter.setPointSets(sampleSet);
		scatter.setPointOpacity(255/3);
		scatter.setBackground(Color.white);
		scatter.display("UAPCA-test");
		
		
		UAPCA<M> uapca = new UAPCA<>(mc);
		M uaCov = uapca.calcUACov(dataset);
		M uaCovV2 = uapca.calcUACovV2(dataset);
		
		ArrayList<M> sampleList = new ArrayList<>();
		int nSamples = 500_000;
		for(int i=0; i< dataset.size(); i++) {
			NRV<M> nrv = dataset.get(i);
			M samp = nrv.drawSamples(nSamples);
			sampleList.add(samp);
		}
		M covSum = null;
		for(int i=0; i<nSamples; i++) {
			int i_=i;
			M batch = sampleList.stream().map(samp -> mc.getRow(samp, i_)).reduce(mc::concatVert).get();
			M cov = NRV.estimateFromData(mc, batch).cov;
			covSum = covSum==null? cov : mc.add(covSum, cov);
		}
		mc.scale_inp(covSum, 1.0/nSamples);
		
		M sampleCov = NRV.estimateFromData(mc, sampleList.stream().reduce(mc::concatVert).get() ).cov;
		 
		
		System.out.println("ua cov");
		System.out.println(uaCov);
		System.out.println("ua cov V2");
		System.out.println(uaCovV2);
		System.out.println("sample cov");
		System.out.println(sampleCov);
		System.out.println("sample cov sum");
		System.out.println(covSum);
		
	}
	
	static <M> M cov2Drot(MatCalc<M> mc, M diagV, double rad) {
		double cos = Math.cos(rad);
		double sin = Math.sin(rad);
		M rot = mc.matOf(2, cos, -sin, sin, cos);
		return mc.mult_abcT(rot, mc.diagM(diagV), rot);
	}
	
}
