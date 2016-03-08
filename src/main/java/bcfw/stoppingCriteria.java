package bcfw;

public class stoppingCriteria{
	public stoppingCriteria(double bRG, int maxIterations) {
		super();
		BRG = bRG;
		this.maxIterations = maxIterations;
	}
	double BRG; // Best Relative Gap , suggestion  = 1.0E-04
	int maxIterations;
	public boolean shouldStop(stoppingCriteria obtained){
		if(obtained==null) return false;//Zeroth step
		return(BRG>obtained.BRG) || (maxIterations<obtained.maxIterations);
	}
}