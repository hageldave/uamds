package datasets;

import java.util.stream.IntStream;
import java.util.stream.Stream;

import optimization.generic.numerics.MatCalc;
import uamds.NRV;

import static uamds.Utils.sq;

public class StudentGrades {
	public static final double var_verybad = sq(6-0)/12;
	public static final double var_bad = sq(7-6)/12;
	public static final double var_fairlybad = sq(10-5)/12;
	public static final double var_fairlygood = sq(14-10)/12;
	public static final double var_good = sq(18-13)/12;
	public static final double var_verygood = sq(20-14)/12;
	
	public static String[] studentNames = {
		"Tom",
		"David",
		"Bob",
		"Jane",
		"Joe",
		"Jack"
	};
	// placeholder values for textual description "very good", "fairly bad" and so on
	static double[] mark = {-1,-2,-3,-4,-5,-6};
	// replacement means and variances for placeholders
	public static double[][] mark2distrib_1 = { 
			{2, sq(2.0/3)},
			{4.5, sq(5)/18},
			{8, sq(5)/12},
			{12, sq(5)/12},
			{15.5, sq(5)/18},
			{18, sq(2.0/3)}
	};
	public static double[][] mark2distrib_2 = { 
			{4, sq(8)/12},
			{7, sq(7)/12},
			{9, sq(6)/14},
			{13, sq(5)/16},
			{16.5, 1.0},
			{18.2, sq(0.6)}
	};
	
	static double[][] means = {
			{15, mark[3], 14, 15},
			{9, mark[4], mark[3], 10},
			{6, 10.5, 16.5, mark[4]},
			{mark[3], mark[5], 19, 11},
			{mark[0], mark[2], 12, 14},
			{1, 5, 9, 7.5},
	};

	static double[][][] covs = {
			{
				{0,0,0,0},
				{0, mark[3],0,0},
				{0,0,33.33333333333333,0},
				{0,0,0,sq(16-14)/12}
			},
			{
				{0,0,0,0},
				{0,mark[4],0,0},
				{0,0,mark[3],0},
				{0,0,0,0.1}
			},
			{
				{0,0,0,0},
				{0,sq(11-10)/12,0,0},
				{0,0,sq(20-13)/12,0},
				{0,0,0,mark[4]}
			},
			{
				{mark[3],0,0,0},
				{0,mark[5],0,0},
				{0,0,0,0},
				{0,0,0,sq(12-10)/12}
			},
			{
				{mark[0],0,0,0},
				{0,mark[2],0,0},
				{0,0,sq(14-10)/12,0},
				{0,0,0,0}
			},
			{
				{0,0,0,0},
				{0,sq(6-4)/12,0,0},
				{0,0,0,0},
				{0,0,0,sq(9-6)/12}
			}
	};
	
	static double[][] getMeans(double[][] distribValues){
		double[][] m = new double[means.length][means[0].length];
		for(int i=0; i<means.length; i++) {
			for(int j=0; j<means[0].length; j++) {
				double v = means[i][j];
				if(v<0) {
					int markIdx = -((int)v)-1;
					v = distribValues[markIdx][0];
				}
				m[i][j] = v;
			}
		}
		return m;
	}
	
	static double[][][] getCovs(double[][] distribValues){
		double[][][] c = new double[covs.length][covs[0].length][covs[0][0].length];
		for(int i=0; i<covs.length; i++) {
			for(int j=0; j<covs[0].length; j++) {
				for(int k=0; k<covs[0][0].length; k++) {
					double v = covs[i][j][k];
					if(v<0) {
						int markIdx = -((int)v)-1;
						v = distribValues[markIdx][1];
					}
					c[i][j][k] = v;
				}
			}
		}
		return c;
	}
	
	@SuppressWarnings("unchecked")
	public static <M> NRV<M>[] get(MatCalc<M> mc, int which) {
		double[][] mark2Distribs = which==0 ? mark2distrib_1 : mark2distrib_2;
		double[][] means = getMeans(mark2Distribs);
		double[][][] covs = getCovs(mark2Distribs);
		Stream<NRV<M>> stream = IntStream.range(0, means.length)
				.mapToObj(i->{
					return new NRV<>(mc, mc.vecOf(means[i]), mc.matOf(covs[i]));
				});
		return stream.toArray(NRV[]::new);
	}
}
