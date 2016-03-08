package dataStructure;
/***
 * time efficient heap , gives back the minimum
 * keys are positive integers , values are double
 * n is the next empty cell
 * @author jose.ferraro
 *left child is 2*i +1
 *right child is 2*i+2
 *father is (i-1)/2
 */
public class Heap {
	int[] keys;
	double[] values;
	int capacity;
	public int n;
	public Heap(int capacity ){
		this.capacity = capacity;
		keys = new int[capacity];
		values = new double[capacity];
		n = 0;
	}
	
	public void put(int k , double v){
		keys[n] = k;
		values[n] = v;
		n++;
		rise(n-1);
	}
	/**
	 * returns the least value key and removes it from the heap 
	 * @return
	 */
	public int pop(){
		if(n==0) return -1;
		int ans = keys[0];
		delete(0);
		return ans;
	}
	/***
	 * deleta o elemento na posição i e limpa o heap
	 * supposes n>=1
	 * @param i smaller than n
	 */
	private void delete(int i) {
		if(i==n-1){//deleting the last element
			n--;
			return;
		}else if(2*i+2>=n){//does no have children or has only left child
			keys[i] = keys[n-1];
			values[i] = values[n-1];
			n--;
			rise(i);
			return;
		}else{//has both childs
			int smaller_child = 2*i+1; //left one by default
			if(values[smaller_child]>values[smaller_child+1]){//the right one is smaller
				smaller_child++;
			}
			keys[i] = keys[smaller_child];
			values[i] = values[smaller_child];
			delete(smaller_child);
		}
	}

	private void rise(int i) {
		if(i==0){
			return;
		}
		int father = (i-1)/2;
		if(values[father]>values[i]){//order is inversed
			double vi = values[i];
			int ki = keys[i];
			values[i] = values[father];
			keys[i] = keys[father];
			values[father] = vi;
			keys[father] = ki;
			rise(father);
		}
	}

	public String str() {
		String ans = "{";
		for(int i = 0 ; i<n ; i++){
			ans = ans+" "+values[i]+" ";
		}
		return ans+"}";
	}
	
}
