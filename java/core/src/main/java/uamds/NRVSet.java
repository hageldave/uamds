package uamds;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Collection of normally distributed random vectors ({@link NRV}s).
 * This is the class used to represent uncertainty-aware data for UAMDS.
 * 
 * @param <M> matrix data type
 */
public class NRVSet<M> extends ArrayList<NRV<M>> {
	private static final long serialVersionUID = 1L;

	public M mean(int i) {
		return get(i).mean;
	}
	
	public M cov(int i) {
		return get(i).cov;
	}
	
	public NRVSet<M> copy() {
		NRVSet<M> c = new NRVSet<>();
		for(int i=0; i<size(); i++) {
			c.add(get(i).copy());
		}
		return c;
	}
	
	public static <M> NRVSet<M> fromArray(NRV<M>[] nrvs) {
		NRVSet<M> set = new NRVSet<>();
		set.addAll(Arrays.asList(nrvs));
		return set;
	}
}
