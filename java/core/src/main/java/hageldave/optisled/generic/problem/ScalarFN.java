package hageldave.optisled.generic.problem;

import hageldave.optisled.generic.numerics.MatCalc;

/**
 * Function taking vector input and giving scalar output.
 * @param <M> matrix (vector) type
 */
public interface ScalarFN<M> {
	/**
	 * @param x function argument (vector)
	 * @return value (scalar) of function evaluated at x
	 */
	public double evaluate(M x);

	/**
	 * Special {@link ScalarFN} with gradient information.
	 * Gradients should be column vectors (by convention), e.g.
	 * {@code d/dx x^T A x = A^T 2x} and NOT {@code 2x^T A}, or using a linear example,
	 * {@code d/dx a^T x = d/dx dot(a,x) = a} and NOT {@code a^T}.
	 * @param <M> matrix type
	 */
	public interface ScalarFNWithGradient<M> extends ScalarFN<M> {
		public VectorFN<M> gradient();
	}
	
	public static <M> ScalarFNWithGradient<M> constant(MatCalc<M> mc, final double c) {
		return new ScalarFNWithGradient<M>() {
			@Override
			public double evaluate(M x) {
				return c;
			}
			@Override
			public VectorFN<M> gradient() {
				return x->mc.scale(x, 0.0);
			}
		};
	}
	
	public static <M> ScalarFNWithGradient<M> linear(final MatCalc<M> mc, final M coefficients, final double c) {
		return new ScalarFNWithGradient<M>() {
			@Override
			public double evaluate(M x) {
				return mc.inner(coefficients,x) + c;
			}
			@Override
			public VectorFN<M> gradient() {
				return x->coefficients;
			}
		};
	}

	public static <M> ScalarFNWithGradient<M> quadratic(MatCalc<M> mc, final M quad, final double c){
		return quadratic(mc,quad,null,c);
	}

	public static <M> ScalarFNWithGradient<M> quadratic(MatCalc<M> mc, final M quad, final M lin, final double c) {
		if (lin != null)
			return new ScalarFNWithGradient<M>() {
				@Override
				public double evaluate(M x) {
					return mc.inner(x, mc.mult_ab(quad, x)) + mc.inner(lin, x) + c;
				}

				@Override
				public VectorFN<M> gradient() {
					M quadT = mc.trp(quad);
					return new VectorFN<M>() {
						@Override
						public M evaluate(M x) {
							return mc.add( mc.scale(mc.mult_ab(quadT,x),2) , lin );
						}
					};
				}
			};
		else
			return new ScalarFNWithGradient<M>() {
				@Override
				public double evaluate(M x) {
					return mc.inner(x, mc.mult_ab(quad, x)) + c;
				}

				@Override
				public VectorFN<M> gradient() {
					return new VectorFN<M>() {
						@Override
						public M evaluate(M x) {
							return mc.scale(mc.mult_ab(quad,x),2);
						}
					};
				}
			};
	}
	
}