package uamds;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class RVPointSet<M> extends ArrayList<NRV<M>> {
	private static final long serialVersionUID = 1L;

	public M mean(int i) {
		return get(i).mean;
	}
	
	public M cov(int i) {
		return get(i).cov;
	}
	
	public RVPointSet<M> translateMean(ArrayList<M> deltaMeans, double m) {
		RVPointSet<M> translated = IntStream.range(0, size())
			.mapToObj(i -> get(i).translateMean(deltaMeans.get(i), m))
			.collect(Collectors.toCollection(RVPointSet::new));
		return translated;
	}
	
	public RVPointSet<M> copy() {
		RVPointSet<M> c = new RVPointSet<>();
		for(int i=0; i<size(); i++) {
			c.add(get(i).copy());
		}
		return c;
	}
	
	public static <M> RVPointSet<M> fromArray(NRV<M>[] nrvs) {
		RVPointSet<M> set = new RVPointSet<>();
		set.addAll(Arrays.asList(nrvs));
		return set;
	}
}
