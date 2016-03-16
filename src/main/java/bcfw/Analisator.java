package bcfw;


/***
 * the class that provides the results of an assignment
 */

public interface Analisator<T> {
	/**
	 *
	 * @param path_matrix is such that path_matrix[a][b][c] is the c-th link on the path from node a to node b
	 * @param od_matrix the origin destination matrix
	 * @return
	 */
	public T simplePath(int[][][] path_matrix , double[][] od_matrix);
	/**
	 * combines two intermediate results into one
	 * @param left
	 * @param right
	 * @param left_proportion
	 * @return
	 */
	public T combine(T left,T right, double left_proportion);
}
