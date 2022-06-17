package optimization.generic.solver;

import optimization.generic.numerics.MatCalc;
import optimization.generic.problem.ScalarFN;
import optimization.generic.problem.VectorFN;
import optimization.history.DescentLog;

/**
 * Gradient descent instance with hyper parameter configuration.
 * Descent is run by calling {@link #arg_min(ScalarFN, VectorFN, M, DescentLog)}.
 * @param <M> matrix type
 */
public class GradientDescent<M> {

	/** 
	 * factor by which the step size is decreased during line search, 
	 * in ]0,1[ 
	 */
	public double stepDecr = 0.5;
	/** 
	 * factor by which the step size is increased again after a step has 
	 * been taken, in [1,infty] 
	 */ 
	public double stepIncr = 1.2;
	/** 
	 * factor for determining sufficient decrease during line search 
	 * (see 1st wolfe condition), in [0.01,0.1] typically 
	 */
	public double lineSearchFactor = 0.01;
	/**
	 * initial step size the algorithm starts with
	 */
	public double initialStepSize = 1.0;
	/**
	 * when the algorithm's steps have decreased below this step size threshold
	 * it terminates, thinking it has reached the minimum
	 */
	public double terminationStepSize = 1e-8;
	/**
	 * maximum number of line search iterations 
	 * (preventing infinite loops in ill conditioned problems) 
	 */
	public int maxLineSearchIter = 100;
	/**
	 * maximum number of descent steps to take
	 * (preventing infinite loops in ill conditioned problems) 
	 */
	public int maxDescentSteps = 100;

	/** the step size of gd when argmin terminated */
	public double stepSizeOnTermination;

	/** the matrix calculation object for the matrix type M */
	public final MatCalc<M> mc;

	/** the loss when arg_min terminates */
	public double lossOnTermination;

	/**
	 * Creates a new GD instance for matrices of type M using
	 * specified matrix calculator.
	 * @param mc matrix calculator to perform linear algebra calculations
	 */
	public GradientDescent(MatCalc<M> mc) {
		this.mc = mc;
	}

	/**
	 * finds argmin by performing gradient descent
	 * @param f function to be minimized
	 * @param df gradient of the function
	 * @param initialGuess initialization (guess of minimum location)
	 * @return location of minimum
	 */
	public M arg_min(ScalarFN<M> f, VectorFN<M> df, M initialGuess){
		return this.arg_min(f,df,initialGuess,null);
	}

	/**
	 * finds argmin by performing gradient descent
	 * @param f function to be minimized
	 * @param df gradient of the function
	 * @param initialGuess initialization (guess of minimum location)
	 * @param log optional log object for recording the optimization trajectory (can be null)
	 * @return location of minimum
	 */
	public M arg_min(ScalarFN<M> f, VectorFN<M> df, M initialGuess, DescentLog log){
		double a = initialStepSize;
		M x = mc.copy(initialGuess);
		int numSteps = 0;
		//
		double fx;
		M dfx;
		M d;
		M step;
		do {
			fx = f.evaluate(x);
			dfx = df.evaluate(x);
			d = mc.normalize_inp(mc.scale(dfx, -1.0));
			if(log != null) {
				log.position(mc.toArray(x));
				log.loss(fx);
				log.direction(mc.toArray(d));
				log.stepSize(a);
			}
			// perform line search
			int numLinsrchIter = 0;
			// while( f(x+a*d) > f(x) + df(x)'a*d*l ) 1st wolfe condition
			while( 
					f.evaluate(mc.add(x, step=mc.scale(d,a))) > fx + mc.inner(dfx,step)*lineSearchFactor
					&& numLinsrchIter++ < maxLineSearchIter
			){
				a *= stepDecr;
				if(log != null)
					log.stepSize(a);
			}
			// update location
			x = mc.add(x,step);
			stepSizeOnTermination = a;
			a *= stepIncr;
		} while( ++numSteps < maxDescentSteps && mc.norm(step) > terminationStepSize );

		this.lossOnTermination = f.evaluate(x);
		if(log != null) {
			log.position(mc.toArray(x));
			log.loss(lossOnTermination);
		}
		
		return x;
	}
	
}
