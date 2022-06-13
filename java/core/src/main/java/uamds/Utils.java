package uamds;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.net.URL;
import java.time.LocalDateTime;
import java.time.temporal.ChronoUnit;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.ForkJoinPool;
import java.util.function.Supplier;
import java.util.stream.IntStream;

import hageldave.optisled.generic.numerics.MatCalc;

public class Utils {

	public static void requireEquals(Object a, Object b, Supplier<String> err) {
		if(a == null) {
			if(a!=b)
				throw new RuntimeException(err.get());
		} else {
			if(!(a == b || a.equals(b))) {
				throw new RuntimeException(err.get());
			}
		}
	}
	
	/** x^2 */
	public static double sq(double x) {
		return x*x;
	}
	
	
	public static String arrayToString(double[] arr, int decimalPlaces) {
		String format = "%."+decimalPlaces+"f";
		StringBuilder sb = new StringBuilder();
		sb.append("[");
		for(int i=0; i<arr.length; i++) {
			sb.append(String.format(Locale.US, format, arr[i]));
			if(i < arr.length-1)
				sb.append(", ");
		}
		sb.append("]");
		return sb.toString();
	}
	
	public static <M> M calcProjectionSubspaceSimilarities(MatCalc<M> mc, M[] projections) {
		M[] normalizedProjs = Arrays.stream(projections)
				.map(p->{
					M lengths = mc.sqrt_inp(mc.rowSums(mc.elmmul(p, p)));
					M divByLenghts = mc.elemwise_inp(lengths, v->1.0/v);
					return  mc.mulRowsByColVec(p, divByLenghts);
				})
				.toArray(mc::matArray);
		double[][] similarities = new double[projections.length][projections.length];
		for(int i=0; i<projections.length; i++) {
			for(int j=i; j<projections.length; j++) {
				similarities[i][j] = similarities[j][i] = calcProjectionSubspaceSimilarity(mc, normalizedProjs[i], normalizedProjs[j]);
			}
		}
		return mc.matOf(similarities);
	}
	
	public static <M> double calcProjectionSubspaceSimilarity(MatCalc<M> mc, M p1, M p2) {
		M dots = mc.mult_abT(p1, p2);
		return mc.norm(dots)/Math.sqrt(mc.numRows(p1));
	}
	
	public static void writeObjToFile(Object obj) {
		String filename = LocalDateTime.now().truncatedTo(ChronoUnit.SECONDS).toString().replace(':', '-')+".obj";
		ForkJoinPool.commonPool().execute(()->{
			try(
					FileOutputStream fos = new FileOutputStream(filename);
					BufferedOutputStream bos = new BufferedOutputStream(fos);
					ObjectOutputStream oos = new ObjectOutputStream(bos);
					){
				oos.writeObject(obj);
			} catch (IOException e) {
				e.printStackTrace();
			}
		});
	}
	
	public static <T> T loadObjFromFile(URL url, Class<T> type) {
		try(
				InputStream is = url.openStream();
				BufferedInputStream bis = new BufferedInputStream(is);
				ObjectInputStream ois = new ObjectInputStream(bis);
		){
			Object obj = ois.readObject();
			return type.cast(obj);
		} catch (IOException | ClassNotFoundException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	public static int[] argsort(final double[] toSort) {
		Integer[] indices = IntStream.range(0, toSort.length).mapToObj(Integer::valueOf).toArray(Integer[]::new);
		Arrays.sort(indices, (i,j)->Double.compare(toSort[i], toSort[j]));
		return Arrays.stream(indices).mapToInt(Integer::intValue).toArray();
	}
	
}
