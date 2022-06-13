package hageldave.optisled.generic.solver;

public class TrajectoryInfo implements Cloneable {
	public double[] x;
	public double fx;
	public double[] gx;
	public double[] lambda;
	public double loss;
	public double mu;
	public boolean isGradientDescent;
	
	public TrajectoryInfo copy() {
		try {
			return (TrajectoryInfo)this.clone();
		} catch (CloneNotSupportedException e) {
			// cannot happen since we are cloneable
			return null;
		}
	}
}
