package uamds.other;

import java.util.LinkedList;
import java.util.function.Consumer;

/** Pointer/Reference class. With notification ability */
public final class Ref<T> {

	public T r;
	
	private LinkedList<Consumer<T>> listeners;
	
	public Ref(T r) {
		this.r=r;
	}
	
	public Ref() {}
	
	public T get() {
		return r;
	}
	
	public void set(T r) {
		this.r = r;
		notifyListeners();
	}
	
	public boolean isNull() {
		return r == null;
	}
	
	public static <T> Ref<T> of(T r){
		return new Ref<>(r);
	}
	
	public synchronized Consumer<T> addListener(Consumer<T> l){
		if(listeners == null)
			listeners = new LinkedList<>();
		listeners.add(l);
		return l;
	}
	
	public synchronized void removeListener(Consumer<T> l) {
		if(listeners != null)
			listeners.remove(l);
	}
	
	private synchronized void notifyListeners() {
		if(listeners != null)
			listeners.forEach(l->l.accept(get()));
	}
}
