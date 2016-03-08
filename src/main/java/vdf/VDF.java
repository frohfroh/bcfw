package vdf;

public abstract class VDF {
	public abstract String[] attributesNames();
	
	public abstract double valueAt(double d, double[] is) ;

	public abstract double derivative(double d, double[] is) ;

	public abstract double integral(double d, double[] is) ;

}
