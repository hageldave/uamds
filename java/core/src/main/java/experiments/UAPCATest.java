package experiments;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.stream.IntStream;

import uamds.datasets.StudentGrades;
import uamds.demo.LoFiScatter;
import hageldave.optisled.ejml.MatCalcEJML;
import hageldave.optisled.generic.numerics.MatCalc;
import uamds.other.NRV;
import uamds.other.NRVSet;
import uamds.other.UAPCA;

public class UAPCATest {

	public static void main(String[] args) {
//		testcov(new MatCalcEJML());
		testDensityCovCalc(new MatCalcEJML());
	}
	
	static <M> void testcov(MatCalc<M> mc) {
		NRV<M> normal4d = NRV.randNRV(mc, 4);
		
		M samples4d = normal4d.drawSamples(10_000_000);
		M samples2d1 = mc.getRange(samples4d, 0, mc.numRows(samples4d), 0, 2);
		M samples2d2 = mc.getRange(samples4d, 0, mc.numRows(samples4d), 2, 4);
		
		NRV<M> nrv2d1 = estimateFromData(mc, samples2d1);
		NRV<M> nrv2d2 = estimateFromData(mc, samples2d2);
		NRV<M> reconst4d = estimateFromData(mc, samples4d);
		
		M block1st = mc.getRange(normal4d.cov, 0, 2, 0, 2);
		
		System.out.println(block1st);
		System.out.println(nrv2d1.cov);
		System.out.println(mc.frob(mc.sub(block1st, nrv2d1.cov)));
		
		System.out.println(normal4d.cov);
		System.out.println(reconst4d.cov);
		System.out.println(mc.frob(mc.sub(reconst4d.cov, normal4d.cov)));
		
		System.out.println(mc.frob(mc.sub(mc.getRange(reconst4d.cov, 0, 2, 0, 2), nrv2d1.cov)));
	}
	
	static <M> void testDensityCovCalc(MatCalc<M> mc) {
		NRV<M> nrv = new NRV<>(mc, 4);
		int n = 1_000;
		M samples = nrv.drawSamples(n);
//		M densities = nrv.pdf().evalManyAt(samples);
		M densities = mc.add(mc.zeros(n), 1.0/n);
		NRV<M> nrvEst1 = estimateFromData(mc, samples);
		NRV<M> nrvEst2 = estimateFromData(mc, samples, densities);
		
		System.out.println(nrvEst1.cov);
		System.out.println(nrvEst2.cov);
		System.out.println(mc.frob(mc.sub(nrvEst1.cov, nrvEst2.cov)));
	}
	
	public static <M> NRV<M> estimateFromData(MatCalc<M> mc, M data) {
		M mean = mc.colMeans(data);
		data = mc.subRowVec(data, mean);
		M cov = mc.matmul(mc.trp(data), data);
		cov = mc.scale_inp(cov, 1.0/mc.numRows(data));
		
		return new NRV<M>(mc, mc.trp(mean), cov);
	}
	
	public static <M> NRV<M> estimateFromData(MatCalc<M> mc, M data, M densities) {
		int rows = mc.numRows(data), cols = mc.numCols(data);
		M cov = mc.zeros(cols, cols);
		
		M mean = mc.colMeans(data);
		data = mc.subRowVec(data, mean);
		data = mc.mulRowsByColVec(data, densities);
		
		for(int i=0; i<rows; i++) {
			M row = mc.getRow(data, i);
			M outerProducts = outer(mc, row, data);
			mc.add_inp(cov, outerProducts);
		}
		return new NRV<M>(mc, mean, cov);
	}
	
	/**
	  * @param a vector (number of elements equals number of columns of b)
	  * @param b matrix/row-vector (number of columns equals number of elements of a)
	  * @return sum of outer products. sum_i a*b_i where b_i is the ith row of b.
	  */
	 static <M> M outer(MatCalc<M>mc, M a, M b) {
		 M colSums = mc.colSums(b);
		 return IntStream.range(0, mc.numElem(a))
		 .mapToObj(i->mc.scale(colSums, mc.get(a, i)))
		 .reduce(mc::concatVert)
		 .get();
	 }
	
	public static <M> M shuffleRows(MatCalc<M> mc, M m) {
		ArrayList<double[]> rows = new ArrayList<>(Arrays.asList(mc.toArray2D(m)));
		Collections.shuffle(rows);
		return mc.matOf(rows.toArray(double[][]::new));
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
