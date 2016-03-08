package bcfw;


/***
 * matriz_de_caminhos[a][b][c] é o c-ésimo link no caminho do nó a ao nó b
 */

public interface Analisator<T> {
	public T simplePath(int[][][] matriz_de_caminhos , double[][] matriz_od);
	public T combine(T left,T right, double left_proportion);
}
