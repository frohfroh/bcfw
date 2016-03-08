package bcfw;

public interface SimpleAnalisator<T> extends FinishingAnalisator<T, T> {
	default T finalize(T bruto){
		return bruto;
	}
	
	
}
