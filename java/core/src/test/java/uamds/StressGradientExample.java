package uamds;

import java.util.stream.IntStream;

import hageldave.optisled.ejml.MatCalcEJML;
import hageldave.optisled.generic.numerics.MatCalc;
import uamds.UAMDS.PreCalculatedValues;
import uamds.other.NRV;
import uamds.other.NRVSet;

public class StressGradientExample {

	
	public static void main(String[] args) {
		test(new MatCalcEJML());
	}
	
	static <M> M mk_cov_mat(MatCalc<M> mc, int d, double s) {
		M a = mc.matOf(d, IntStream.range(0, d*d).mapToDouble(i-> Math.sqrt(i*s)).toArray());
		a = mc.subRowVec(a, mc.colMeans(a));
		M cov = mc.mult_aTb(a, a);
		return mc.scale(cov, 1.0/d);
		//return mc.eye(d, s);
	}
	
	static <M> void test(MatCalc<M> mc) {
		int n = 4;
		int d = 6;
		int lo_d = 2;
		
		NRVSet<M> set = new NRVSet<>();
		M[] c = mc.matArray(n);
		M[] B = mc.matArray(n);
		for(int i=0; i<n; i++) {
			set.add(new NRV<M>(mc, mc.scale(mc.ones(d), i+1), mk_cov_mat(mc, d, i+1)));
			c[i] = mc.scale(mc.ones(lo_d), i+1);
			B[i] = mc.scale(mc.ones(lo_d, d), i*0.1+0.1);
		}
		
		PreCalculatedValues<M> pre = new PreCalculatedValues<>(mc, set);
		
		double s = UAMDS.stressFromProjecton_ij(mc, 1, 2, pre, B, c, d, lo_d, null, null);
		M[][] grad = UAMDS.gradientFromProjection_ij(mc, 1, 2, pre, B, c, d, lo_d);
		
		System.out.println(s);
		System.out.println(UAMDS.stressFromProjecton(mc, pre, B, c, d, lo_d, null, null));
		
		
		
	}
	
}
