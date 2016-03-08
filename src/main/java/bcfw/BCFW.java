package bcfw;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.PriorityQueue;
import java.util.Random;




//import teste.WPMainfunction;
import vdf.VDF;
public class BCFW {
	
	public BCFW(stoppingCriteria sc){
		this.sc = sc;
	}
	// O que segue é um grafo
	public int V; //número de nós
	public int E;//número de links
	public int Centroides; //os nós de 0 a Centroides-1 são centróides
	public int fromNode[]; //node form de cada ligação
	public int toNode[]; //node To dde cada ligação
	public int arrivingLinks[][];//ligações entrantes
	public int exitingLinks[][];//ligações saintes
	public double linkParameters[][];//parâmetros de cada ligação
	public  VDF[] VDFs;	//Funções vdfs
	public boolean useConnectors = true;
	long p1 , p2,p3,p4;
	//O que segue é a matriz 
	public double[][] OD;
	
	//Constantes
	final int LINE_SEARCH_POINTS = 12;//number of points to evaluate
	final double LINE_SEARCH_GROWTH = 1.2;// how much can the step grow each step
	final stoppingCriteria sc;
	
	class step_result{
		double β_0;
		double β_1;
		double β_2;
		double τ;
		double[] d;//It is not necessary as one could recalculate, but it can save time(We hope), no idea if the compiler could optimize it
		int[][][] caminhos;
		double lowerBound;
		public step_result(double β_0, double β_1, double β_2, double τ, double[] d, int[][][] caminhos, double lowerBound) {
			super();
			this.β_0 = β_0;
			this.β_1 = β_1;
			this.β_2 = β_2;
			this.τ = τ;
			this.caminhos = caminhos;
			this.d = d;
			this.lowerBound = lowerBound;
		}
	}
	
	public <T> T assign(Analisator<T> analisator){
		int k = 1;
		double tot=0.0;
		for(double[] arr : OD) for (double d : arr) tot+=d;
		System.out.println("A matriz OD JAVA original tem "+tot);
		double[] empty_network = new double[E];//Full of zeros by default 
		int[][][] shortest_paths =  ShortestPath(empty_network);
		double[] x_k = putVolumes(shortest_paths);
		//The variables ending in "r" are the result asked
		T s_km2r = analisator.simplePath(shortest_paths, OD);
		T x_kr = s_km2r;
		T s_km1r = s_km2r , s_kr = s_km2r ;//initializes the results
		T x_kp1r=null;
		double[] s_km1 = null,s_km2 = null;
		double τ_km1 = 1,τ_km2 = 1;
		double BLB = Double.NEGATIVE_INFINITY;//Best lower bound
		stoppingCriteria obtained = null;
		step_result r = null;
		p1 = 0;
		p2 = 0;	
		p3 = 0;
		p4 = 0;

		while(!sc.shouldStop(obtained)){
			System.out.println("Iteração "+k+" objetivo: "+ObjectiveFunction(x_k));//double calculation of ObjectiveFunction
			/*System.out.println("Les volumes sont: ");
			for(int i =0; i<x_k.length;i++) System.out.print(" "+x_k[i]);*/
			//System.out.println("novos volumes");
			//System.out.println("2392 "+x_k[2392]);
			//System.out.println("712 "+x_k[712]);
			//System.out.println("2433 "+x_k[2433]);

			System.out.print("\n");
			r = normal_step(x_k,s_km1,s_km2,τ_km1,τ_km2);

			double[] d_k = r.d;
			double τ_k = r.τ;
			double[] x_kp1 = Positive(Plus(x_k ,Times(τ_k  , d_k)));//Should always be positive except for numeric errors
			double[] s_k = Plus(x_k , d_k);

			System.out.println("o passo é "+τ_k);
			/*System.out.println("A direção dos novos volumes é ");
			for(double vol : d_k) System.out.print(vol+" ");*/
			System.out.println("");		

			//This formula explains the analysis
			//β_0+β_1+β_2 = 1
			//s_k = β_0 * y_k + β_1*s_km1 +β_2*s_km2 = β_0 * y_k + (1-β_0)*(β_1/(β_1+β_2)*s_km1 + β_2/(β_1+β_2)*s_km2)

			
			
			T y_kr = analisator.simplePath(r.caminhos , OD);

			if(r.β_1+r.β_2>0){
				T right = analisator.combine(s_km1r , s_km2r,r.β_1/(r.β_1+r.β_2) );
				s_kr = analisator.combine(y_kr, right,  r.β_0);
			}else{//if we have a simple FW step
				s_kr = y_kr;
			}
			//x_kp1 = x_k + d_k * τ_k
			//x_kp1 = x_k * (1 - τ_k) + s_k * τ_k
			//x_kp1 = x_k*(1-τ_k) + y_k*β_0*τ_k + s_km1*β_1*τ_k +  s_km2*β_2*τ_k
			x_kp1r = analisator.combine(x_kr , s_kr , (1-τ_k) );


			if(BLB<r.lowerBound) BLB = r.lowerBound;
			//We go to the next iteration
			k++;
			τ_km2 = τ_km1;
			τ_km1 = τ_k;
			x_k = x_kp1;
			s_km2 = s_km1;
			s_km1 = s_k;
			x_kr = x_kp1r;
			s_km2r = s_km1r;
			s_km1r = s_kr;
			double UB = this.ObjectiveFunction(x_k);//Upper Bound
			double BRG = (UB - BLB)/UB;//Best relative Gap
			obtained = new stoppingCriteria(BRG , k);
		}
		System.out.println("\n\n\n Tempos ");
		System.out.println("p1 "+p1/1000000);
		System.out.println("p2 "+p2/1000000);
		System.out.println("p3 "+p3/1000000);
		System.out.println("p4 "+p4/1000000+"\n\n\n");

		//for(int i = 0 ; i< x_k.length;i++) System.out.println(x_k[i]);
		return x_kr;
	}
	//gives x if x>0 0 otherwise
	//used to prevent precision errors, as volumes should be non-negative
	private static double[] Positive(double[] plus) {
		plus = plus.clone();
		for(int i =0; i<plus.length ; i++) if(plus[i]<0.0) plus[i]=0.0;
		return plus;
	}
	//TODO passo negativo
	//TODO passo nulo -> todo o resto é zero
	/**
	 * Os resultados são constantes que devem permitir a um analisador fazer seu trabalho.
	 * No entanto, os novos volumes devem forçosamente ser calculados.
	 * Esta etapa funciona como um Frank-Wolfe normal se os dois passos anteriores não forem estritamente maiores que zero
	 * @param x_k volumes nas ligação da iteração anterior(k)
	 * @return um objeto com os fatores, para se analisar
	 */
	public step_result normal_step(double[] x_k , double[] s_km1 , double[] s_km2,double τ_km1,double τ_km2 ){

		int[][][] matriz_de_caminhos = ShortestPath(x_k);//of a "full-of-ones matrix"

		double[] y_k = putVolumes(matriz_de_caminhos);
		System.out.println("velho tempo "+Gradient(x_k,x_k));
		System.out.println("desejado tempo "+Gradient(x_k,y_k));
		//double τ_k_num = Gradient(x_k,d_k);
		//System.out.println("A derivada anal é "+τ_k_num);

		double[] d_k_FW = Minus(y_k,x_k);

		if(τ_km1>=1.0 || τ_km2>=1.0 || τ_km1<=0.0 || τ_km2<=0.0){ //the negative cases are likely when the problems dimension is so small that by conjugating you get to an ascending direction 
			System.out.println("Linear search will be performed");
			//Scanner scan = new Scanner(System.in);
			//scan.next();
			double τ_k=LineSearch(x_k , d_k_FW,Math.abs(τ_km1)); //Si négatif, on essaye juste d'en garder l'ordre de grandeur
			return new step_result(1,0,0,τ_k,d_k_FW,matriz_de_caminhos,Double.NEGATIVE_INFINITY);//The -∞ is so that we do not converge after a such step
		}

		double[] d_km1_bar = Minus(s_km1 , x_k);
		double[] d_km2_bar_bar = Plus(Minus(Times(τ_km1,s_km1),x_k),Times(1-τ_km1,s_km2));
		double μ_k_num = Hessian(d_km2_bar_bar , x_k , d_k_FW);//TODO save a derivative calculation by deriving first and then calling Hessian?
		double μ_k_den =  Hessian(d_km2_bar_bar , x_k , Minus(s_km2 ,s_km1));
		double μ_k = -μ_k_num/μ_k_den;
		double ν_k_left_num = Hessian(d_km1_bar, x_k ,d_k_FW) ;//Isso(ν) é infelizmente um nu!
		double ν_k_left_den = Hessian(d_km1_bar, x_k ,d_km1_bar) ;
		double ν_k = -ν_k_left_num/ν_k_left_den + μ_k * τ_km1/(1-τ_km1);
		if(μ_k<0.0) μ_k = 0.0; //We must have a convex combination in order to remain feasible
		if(ν_k<0.0) ν_k = 0.0;		
		double β_0_k = 1/(1+ μ_k+ν_k);
		double β_1_k = β_0_k * ν_k;
		double β_2_k = β_0_k * μ_k;
		System.out.println("μ "+μ_k+" ν "+ν_k);
		System.out.println("β_0 "+β_0_k+" β_1 "+β_1_k+" β_2 "+β_2_k);
		double[] d_k = Plus(Plus(Times(β_0_k,d_k_FW),Times(β_1_k,Minus(s_km1,x_k))),Times(β_2_k,Minus(s_km2,x_k)));
		double τ_k_num = Gradient(x_k,d_k);
		double τ_k_den = Hessian(d_k,x_k,d_k);
		double τ_k =  -τ_k_num/τ_k_den;//Newton step
		System.out.println("A derivada anal é "+τ_k_num);

		/*
		 * in this case, the direction odf descent is positive, that is, the objective increases
		 * This happens because when conjugating by the last two directions, it's not given that we get a negative derivative
		 * This(having a negative derivative) should happen most of the time because we just found the minimum along a line in this direction, so the derivative 
		 * in relation to these two directions should be zero, but as the step determination is not perfect , this is not valid
		 * For the direction of iteration n-2 this is even worst, as we're no longer in the same point, so the 
		 * derivative has also changed
		 * We decide to take a FW step, but one coould also optionally choose to try a simple conjugated FW step
		 * This can be tested in  the future
		 */
		if(τ_k < 0.0){
			System.out.println("Linear search will be performed for τ_k < 0.0");
			System.out.println("A derivada anal FW é "+Gradient(x_k,d_k_FW));
			if(Gradient(x_k,d_k_FW) <  0 && Gradient(x_k,d_k_FW) >  0){
				System.out.println("A derivada FW é positiva,  o que não deveria ocorrer nunca");
				try {
					PrintWriter timesF = new PrintWriter("D:\\lixo\\times", "UTF-8");
					double[] times = new double[E];
					for(int i = 0 ; i<E ; i++) times[i] = VDFs[i].valueAt(x_k[i],linkParameters[i]);
					for(int i = 0 ; i<E ; i++) timesF.println(times[i]);
					timesF.close();
					PrintWriter FromF = new PrintWriter("D:\\lixo\\From", "UTF-8");
					for(int i = 0 ; i<E ; i++) FromF.println(this.fromNode[i]);
					FromF.close();
					PrintWriter ToF = new PrintWriter("D:\\lixo\\To", "UTF-8");
					for(int i = 0 ; i<E ; i++) ToF.println(this.toNode[i]);
					ToF.close();
					for(int f = 0 ; f<Centroides ; f++){
						for(int t = 0 ; t<Centroides ; t++){
							PrintWriter pathF = new PrintWriter("D:\\lixo\\path\\"+f+"_"+t, "UTF-8");
							for(int l : matriz_de_caminhos[f][t]){
								pathF.println(l+"\t"+fromNode[l]+"\t"+toNode[l]);
							}
							pathF.close();
					}
					}
					PrintWriter custosF = new PrintWriter("D:\\lixo\\custo", "UTF-8");
					for(int f = 0 ; f<Centroides ; f++){
						for(int t = 0 ; t<Centroides ; t++){
							double cost = 0.0;
							for(int l : matriz_de_caminhos[f][t]) cost+= times[l];
							custosF.println(f+"\t"+t+"\t"+cost+"\t"+OD[f][t]);
					}
					}
					for(int l : matriz_de_caminhos[87][107]){
						System.out.println(l+" "+fromNode[l]+" "+toNode[l]+" "+times[l]);
					}
					custosF.close();

				} catch (FileNotFoundException | UnsupportedEncodingException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			τ_k=LineSearch(x_k , d_k_FW,Math.abs(τ_km1)); //Si négatif, on essaye juste d'en garder l'ordre de grandeur
			
			
			
			
			return new step_result(1,0,0,τ_k,d_k_FW,matriz_de_caminhos,Double.NEGATIVE_INFINITY);//The -∞ is so that we do not converge after a such step

		}
		double τ_k_alt = LineSearch(x_k , d_k,Math.abs(τ_km1));
		double obj_alt = ObjectiveFunction(Plus(x_k , Times(τ_k_alt , d_k)));
		System.out.println("o passo busca linear é "+τ_k_alt+" com objetivo "+obj_alt);

		if(τ_k>1.0) τ_k=1.0;
		double obj = ObjectiveFunction(Plus(x_k , Times(τ_k , d_k)));
		System.out.println("o passo newton é "+τ_k+" com objetivo "+obj);

		if(obj > obj_alt){
			τ_k = τ_k_alt;
		}
		double lowerBound = ObjectiveFunction(x_k)+Gradient(x_k,d_k_FW);
		System.out.println("lower bound "+lowerBound);
		if(τ_k == 0.0){
			double[] ones = new double[E];
			for(int i = 0  ; i<E ; i++) ones[i]=1.0;
			double derivada_dir_d_k = τ_k_num;
			double h = 0.00000000001;
			double numeric = (ObjectiveFunction(Plus(x_k,Times(h ,d_k )))-ObjectiveFunction(x_k))/h;
			System.out.println("analitica "+derivada_dir_d_k+" numérica "+numeric);
			System.out.println("ofbase + \" \"+ofmais+\" \"+numeri+\" \"+anal");
			for(int i = 0 ; i< E ; i++){
				double ofbase =  VDFs[i].integral(x_k[i],linkParameters[i]);
				double ofmais = VDFs[i].integral(x_k[i]+h*d_k[i],linkParameters[i]);
				double numeri =  (ofmais - ofbase)/h;
				double anal = VDFs[i].valueAt(x_k[i],linkParameters[i]) * d_k[i];
				System.out.println(ofbase + " "+ofmais+" "+numeri+" "+anal);
				if(Math.abs(numeri - anal) > 100.0){
					System.out.println(VDFs[i].getClass().getName()+" abo "+VDFs[i].getClass().getSimpleName());
					/*teste.WPMainfunction fun = (WPMainfunction) VDFs[i];*/
					System.out.println(VDFs[i]);
					System.out.println("x x+h " + x_k[i]+"  "+(x_k[i]+h*d_k[i]));
					for( double p : linkParameters[i]) System.out.print(p+ "  ");
					System.out.println("param ");
				}
			}
		}
		return new step_result(β_0_k,β_1_k,β_2_k,τ_k,d_k,matriz_de_caminhos,lowerBound);

	}
	/**
	 * Busca linear para minimizar a função partinfo de x e indo em direção a d.
	 * só é usada nas primeiras etapas ou em casos excepcionais, senão usamos um passo de Newton
	 * @param x
	 * @param d
	 * @param τ_km1 
	 * @return
	 */
	private double LineSearch(double[] x, double[] d, double τ_km1) {
		int lsp = LINE_SEARCH_POINTS;
		double[] searchPoints = new double[lsp];
		double dif = LINE_SEARCH_GROWTH * τ_km1/(lsp-1);
		if(dif*(lsp-1) >1) dif = 1.0/(lsp-1) ;
		for(int i = 0; i< lsp ; i++) searchPoints[i]= i*dif;
		double[] value = new double[lsp];
		for(int i = 0; i< lsp ; i++) {
			double[] volumes = Positive(Plus(x , Times(searchPoints[i],d)));
			value[i] = ObjectiveFunction(volumes);
		}
		System.out.println("The objective and points are: ");
		for(int i =0; i<value.length;i++) System.out.print("("+value[i]+","+searchPoints[i]+") ");
		boolean erro = false;
		
		for(int i =0; i<value.length;i++) if(Double.isNaN(value[i]))erro = true;
		if(erro) for(int i1 = 0; i1<x.length;i1++) System.out.println("x "+x[i1]+" d "+d[i1]);
		
		System.out.println("");		
		double min =  FindMin(value , searchPoints);
		System.out.println("o mínimo é "+min);
		if(min==0.0) {
			if(dif<1.0E-10) return 0.0;
			return LineSearch(x,d,dif);//TODO take this away when we have proper support for Best relative gap stopping criteria
		}
		//we try to increase the step, but if it is already 1.0 there is no point in that
		if(min >= 0.999999 * searchPoints[lsp-1] && min <=0.999999) return LineSearch( x, d, τ_km1*2);
		return min;
	}
	private static double FindMin(double[] value ,double[] searchPoints) {
		double min = Double.POSITIVE_INFINITY;
		int pos = 0;
		for(int i = 0 ; i<value.length;i++){
			if(value[i]<min){
				pos=i;
				min = value[i];
			}
		}
		if(pos==0 || pos==value.length-1) return searchPoints[pos]; //If are extremes, return extreme
		double[] x = {searchPoints[pos-1],searchPoints[pos],searchPoints[pos+1]};
		double[] y = {value[pos-1],value[pos],value[pos+1]};
		return ParabolaMinimum(x,y);
	}
	/**
	 * Givem three points, return the x of the minimum of the parabola
	 * @param x
	 * @param y
	 * @return
	 */
	private static double ParabolaMinimum(double[] x, double[] y) {
		double x1 = x[0];
		double x2 = x[1];
		double x3 = x[2];
		double y1 = y[0];
		double y2 = y[1];
		double y3 = y[2];
		double denom = (x1 - x2)*(x1 - x3)*(x2 - x3);
		double A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
		double B = (x3*x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1*x1* (y2 - y3)) / denom;
		return -B/(2*A);
	}

	public double ObjectiveFunction(double[] volumes) {
		double[] items = new double[E];
		for(int i = 0 ; i<E;i++){
			items[i] = VDFs[i].integral(volumes[i],linkParameters[i]);
		}
		return Sum(items);
	}

	private static double Sum(double[] values) {
		double simples = 0.0;
		for(double v : values) simples+=v;
		double kahan = 0.0;
		double comp = 0.0;
		for(double v : values){
			double y = v -comp;
			double t = kahan + v;
			comp = (t - kahan)-y;
			kahan = t;
		}
		double error = (simples-kahan)/kahan;
		if(kahan==0) error = 0;
		if(error>0) System.out.println("Erro da soma ingênua é "+error);
		return kahan;
	}

	/**
	 * Calcula (∇f)b onde ∇f é o gradiente da função no ponto x
	 * @param x ponto em que se toma o gradiente
	 * @param b vetor a ser multiplicado à direita
	 * @return resultado do produto(vetor)
	 */
	private double Gradient(double[] x, double[] b) {
		double ans = 0;
		for(int i =0; i< x.length;i++){
			double grad = VDFs[i].valueAt(x[i],linkParameters[i]);
			ans+= grad * b[i];
		}
		return ans;
	}
	/** 
	 * Calcula o produto u'Hv , onde H é a hessiana tomada no ponto x
	 * @param u vetor multiplicado à esquerda
	 * @param x ponto em que se toma a Hessiana
	 * @param v vetor multiplicado à direita
	 * @return reultado do produto(real)
	 */
	private double Hessian(double[] u, double[] x, double[] v) {
		double ans = 0.0;
		double[] hessianDiagonal = new double[E];
		for(int i = 0 ; i<v.length ; i++) hessianDiagonal[i] = VDFs[i].derivative(x[i],linkParameters[i]);
		for(int i = 0 ; i<v.length ; i++) ans+=u[i]*hessianDiagonal[i]*v[i];
		return ans;
	}
	private static double[] Plus(double[] u, double[] v) {
		double[] ans = new double[v.length];
		for(int i = 0 ; i<v.length ; i++) ans[i]=u[i]+v[i];
		return ans;
	}
	private static double[] Times(double c, double[] v) {
		double[] ans = new double[v.length];
		for(int i = 0 ; i<v.length ; i++) ans[i]=v[i]*c;
		return ans;
	}
	private static double[] Minus(double[] u, double[] v) {
		double[] ans = new double[v.length];
		for(int i = 0 ; i<v.length ; i++) ans[i]=u[i]-v[i];
		return ans;
	}
	private double[] putVolumes(int[][][] matriz_de_caminhos) {
		double[] volumes = new double[E];
		for(int o = 0; o < Centroides ; o++){
			for( int d = 0 ; d < Centroides ; d++){
				int[] links = matriz_de_caminhos[o][d];
				for(int link : links)	volumes[link]+= OD[o][d];
			}
		}
		return volumes;
	}
	/**
	 * 
	 * @param x volumes na rede
	 * @return rank three "matrix" such that the index are Centroid Centroid numer, 
	 * and the element i j k is the number of the k-th link in the path from Centroid i to Centroid j 
	 */
	private int[][][] ShortestPath(double[] x) {
	
		int[][][] ans = new int[Centroides][][];
		double[] times = new double[E];
		for(int i = 0 ; i<E ; i++) times[i] = VDFs[i].valueAt(x[i],linkParameters[i]);
	
		for(int i = 0 ; i<Centroides ; i++){
			double[] timesConnect = times;
			if (! useConnectors){
				timesConnect = new double[E];
				for(int ii = 0 ; ii<E ; ii++) timesConnect[ii] = times[ii];
				double soma = Sum(times);
				for(int j = 0; j<Centroides ; j++){
					if(j!=i){
						for(int link : exitingLinks[j]){
							timesConnect[link] = soma;
						}
					}
				}
			}
			ans[i] = Dijkstra2(i,timesConnect);//This part can probably be parallelised to huge gains
			
		}

		return ans;
	}
	
	//this is equivalent to Dijkstra, the question is which one is faster
	//this one is also wrong, so before using it would be better to debug
	public int[][] Dijkstra(int root , double[] times) {
		class distancedNode implements Comparable<distancedNode>{
			public distancedNode(int number, double distance , int last , int father) {
				super();
				this.number = number;
				this.distance = distance;
				this.fromLink = last;
				this.fromNode = father;
			}
			int number;
			double distance;
			int fromLink;
			int fromNode;
		@Override
		public int compareTo(distancedNode specified) {
			return (int) Math.signum(this.distance - specified.distance);
		}
		public String toString(){
			return number+" "+distance+" "+fromLink+" "+fromNode;
		}
		}
		double[] distances = new double[E]; //Não deveria ser V?
		for(int n = 0; n<distances.length;n++) distances[n]=Double.POSITIVE_INFINITY;
		distancedNode[] solved_nodes = new distancedNode[V];
		PriorityQueue<distancedNode> candidates = new PriorityQueue<distancedNode>();
		distances[root]=0.0;
		int last_added = root;
		solved_nodes[root] = new distancedNode(root, 0.0 , -1 , root);//he is his own father
		boolean finished = false;
		do{
			double lad = distances[last_added];
			for(int  link : exitingLinks[last_added]){
				double time = times[link];
				int candidate = toNode[link];
				double present_distance = distances[candidate];
				//Nesta linha entrariam os caminhos duplos
				if(lad+time<present_distance)	candidates.add(new distancedNode(candidate,lad+time,link,last_added));
				distances[candidate]=lad+time;
			}
			while(true){
				if(candidates.isEmpty()){
					finished = true;
					break;
				}
				distancedNode next = candidates.poll();
				if(next == null) break;
				if(next.distance> distances[next.number]) continue;//old distancedNode, there is a shorter on candidates
				//following is needed because we can have two paths with equal length to the same node
				if(solved_nodes[next.number] != null) continue;//trying to come back to already visited node
				last_added = next.number;
				solved_nodes[next.number] = next;
				break;
			}
		} while (!finished);
		int paths[][] = new int[Centroides][];
		for(int node =0 ; node<Centroides; node++)	{
			ArrayList<Integer> links = new ArrayList<>();
			distancedNode intermediary_node = solved_nodes[node];
			if(intermediary_node==null) System.out.println("nó é "+node);
			while(intermediary_node.fromLink>=0){
				links.add(intermediary_node.fromLink);
				intermediary_node = solved_nodes[intermediary_node.fromNode];
			}
			int l = links.size();
			int[] path = new int[l];
			for(int i1 = 1 ; i1<=l ; i1++)	path[i1-1] = links.get(l-i1);
			paths[node] = path;
		}

		return paths;
	}
	
	public int[][] Dijkstra2(int root , double[] times){

		int[] tree = DijkstraTree(root ,  times);

		int[][] paths = new int[Centroides][];
		for(int i = 0 ; i  < Centroides ; i++){
			paths[i] = nodes2Links(tree ,i);
		}

		return paths;
	}
	/**
	 * 
	 * @param tree each node point to the node just before in its path
	 * @param i the ending node (the paths ends here)
	 * @return
	 */
	private int[] nodes2Links(int[] tree, int i) {
		int len = tree[2*V+i];//número de links
		int[] nós2 = new int[len+1];
		if(len==0){
			return new int[0];
		}
		int w = 0;
		while(tree[i]!= i){
			nós2[len-w-1] = i; 
			i = tree[i];
			w++;
		}
		nós2[len] = i;

		int[] ans = new int[len];
		for(int j1 =0 ; j1< len ; j1++){
			int tn = nós2[j1];
			ans[j1] = tree[V+tn];
		}
		return ans;
	
	}
	/**
	 * faster Dijkstra
	 * root is whece the paths go
	 * the answer is a  vector of 2 * V sized shortest path tree 
	 * The first V nodes 0...V-1 with the following conventions:
	 * thet root points to itself
	 * unacccessible nodes point to -1
	 * The other nodes point to their father
	 * The nodes V ... 2*V  - 1 are the link that links to its father
	 * The nodes 2*V ... 3*V-1 are the number of links separating it from the root
	 * @param root
	 * @param times
	 * @return
	 */
	
	public int[] DijkstraTree(int root , double[] times) {
		double[] dists = new double[V];
		boolean[] feitos = new boolean[V];
		int[] pais = new int[3*V];
		int number_done = 0;//root is "done" but we still need to add its neighboors 
		pais[root] = root;//the root points to itself
		pais[root+2*V] = 0;
		double inf = Sum(times);
		for(int n = 0; n<dists.length;n++) dists[n]=inf;//uses the sum of all costs as maximum distance
		for(int n = 0; n<dists.length;n++) feitos[n]=false;//uses the sum of all costs as maximum distance
		dists[root] = 0.0;
		dataStructure.Heap heap = new dataStructure.Heap(E); 
		heap.put(root, 0);
		while(number_done < V){
			int next = heap.pop();
			if(next < 0){//no more nodes in the head
				break;//this happens if there are nodes that are not accessible 
			}
			if(feitos[next]) continue;
			feitos[next] = true;
			number_done++;
			for(int link : exitingLinks[next]){
				double parIci = dists[next]+times[link];
				int neighboor = toNode[link];
				if(dists[neighboor]<=parIci) continue;
				if(feitos[neighboor]) continue;
				dists[neighboor] = parIci;
				heap.put(neighboor, parIci);
				pais[neighboor] = next;
				pais[V+neighboor] = link;
				pais[2*V+neighboor] = pais[2*V+next]+1;
			}
		}
		for(int i =0 ; i < V ; i++){//marcamos os nós inacessíveis
			if(!feitos[i]) pais[i] = -1;
		}
		return pais;
	}

}
