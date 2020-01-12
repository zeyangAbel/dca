/******************************************************************************
* This is a class for storing the solution to an instance of DRAP-NC 
****************************************************************************/
public class ResultTypeRAPNC {
	boolean feasible;
	long[] sol;
	public ResultTypeRAPNC(boolean feasible, long[] sol) {
		this.feasible = feasible;
		this.sol = sol;
	}
}
