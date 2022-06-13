package hageldave.optisled.generic.problem;

import hageldave.optisled.generic.numerics.MatCalc;

/**
 * Function taking vector input and giving vector output.
 * @param <M> matrix (vector) type
 */
public interface VectorFN<M> {

	/**
	 * @param x function argument (vector)
	 * @return value (vector) of function evaluated at x
	 */
	public M evaluate(M x);
	
	public static <M> VectorFN<M> constant(final M c){
		return x->c;
	}
	
	public static <M> VectorFN<M> linear(MatCalc<M> mc, final M transform, final M c){
		return x->mc.add(mc.matmul(transform,x),c);
	}
}