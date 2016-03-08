package bcfw;

public interface FinishingAnalisator<T,F> extends Analisator<T>{
	public T simplePath(int[][][] matriz_de_caminhos , double[][] matriz_od);
	public T combine(T left,T right, double left_proportion);
	public F finalize(T bruto);
}
