package uamds.optimization.history;

public interface DescentLog {

	public void position(double[] array);

	public void direction(double[] array);

	public void loss(double fx);

	public void stepSize(double a);

}
